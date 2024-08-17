using CPLEX, JuMP, DelimitedFiles, Random

include("data.jl")
include("SparseMaxFlowMinCut.jl") # lib to solve mincut

ed(i, j) = i < j ? (i, j) : (j, i)

# Distance function for the packing sets
function dist_pack(data, a)
   e = (a[1],a[2])
   if a[1] > a[2]
      e = (a[2],a[1])
   end
   if e[1] == e[2]
      return 0.0
   else
      return data.G′.cost[e]
   end
end


#Classical mctp
function build_model(data::DataCVRP, app::Dict{String,Any})
   pack = Dict()
   E, C, O, M = data.G′.E2, data.C, data.O, data.M
   V = [i-1 for i=1:1 + length(O) + length(M)]

   # Formulation
   mctp = VrpModel()
   # xₑ counts the number of times edge e is used.
   @variable(mctp.formulation, x[e in E], Int)

   # minimizes the total cost
   @objective(mctp.formulation, Min, sum(c(data,e) * x[e] for e in E))

   # ensures that the mandatory facilities are visited exactly once
   @constraint(mctp.formulation, deg1[i in M], sum(x[e] for e in δ′(data, Set([i]))) == 2.0)

   # ensures that the optional facilities are visited at most once
   @constraint(mctp.formulation, deg2[i in O], sum(x[e] for e in δ′(data, Set([i]))) <= 2.0)

   # guarantees that all customers are covered
   @constraint(mctp.formulation, degC[c in C], sum(x[e] for e in δ(data,cover(data,c))) >= 2*data.cov[c])


   # Resource Constrained Shortest Path Graph
   # VRPSolver uses this graph to solve the subproblems by a bi-directional labeling algorithm
   function build_graph()
      V = [i-1 for i=1:1 + length(O) + length(M)]

      # source=sink=depot
      v_source = v_sink = 0
      # minimum and maximum number of paths
      L, U = 0, length(V)-1 
   
      # node ids of G from 0 to n
      G = VrpGraph(mctp, V, v_source, v_sink, (L, U))

      # Resource that corresponds to the capacity of the vehicles
      cap_res_id = add_resource!(G, main = true) # R = R_M = {cap_res_id}
      for i in V
         # Bounds of each vertex
         l_i, u_i = 0.0, Float64(data.p) 
         set_resource_bounds!(G, i, cap_res_id, l_i, u_i)
      end

      # If is length constrained, create another resource for the length
      if app["maxL"] < 10000000.0
         len_res_id = add_resource!(G, main = true)
         max = 0.0
         for i in vcat(O,M)
            e = (0,i)
            if c(data,e) > max
               max = c(data,e)
            end
         end
         for i in V
            l_i, u_i = 0.0, app["maxL"] + 2*max
            set_resource_bounds!(G, i, len_res_id, l_i, u_i)
         end
      end

      E_new = []
      for i in V, j in V
         if i < j
            push!(E_new,(i,j))
         end
      end
      # Creating the arcs of the graph
      for (i,j) in E_new
         e = (i,j)
         
         # arc from i to j
         arc_id = add_arc!(G, i, j)
         # Map the corresponding variable
         add_arc_var_mapping!(G, arc_id, [x[e]])
         if app["maxL"] < 10000000.0
            set_arc_consumption!(G, arc_id, len_res_id, c(data,e))
         end
         if i==0
            set_arc_consumption!(G, arc_id, cap_res_id, 0.5)
         else
            set_arc_consumption!(G, arc_id, cap_res_id, 1.0)
         end

         #arc from j to i
         arc_id = add_arc!(G, j, i)
         add_arc_var_mapping!(G, arc_id, [x[e]])

         if app["maxL"] < 10000000.0
            set_arc_consumption!(G, arc_id, len_res_id, c(data,e))
         end

         if i != 0
            set_arc_consumption!(G, arc_id, cap_res_id, 1.0)
         else
            set_arc_consumption!(G, arc_id, cap_res_id, 0.5)
         end
      end
      return G
   end

   G = build_graph()
   add_graph!(mctp, G)

   # Optional and mandatory nodes are packing sets
   # Packing set is a set that can be visited at most once at an optimal solution
   pack = [v for v in V if v != 0]
   set_vertex_packing_sets!(mctp, [[(G,i)] for i in pack])

   # Define the distance between packing sets (useful for the ng-path relaxation)
   define_elementarity_sets_distance_matrix!(mctp, G, [[dist_pack(data, (i, j)) for j in pack] for i in pack])
   
   # Define the branching variables
   set_branching_priority!(mctp, "x", 1)
   
   # Section for the robust cuts
   cov = [(i, j) for (i, j) in keys(data.G′.cov_matrix) if data.G′.cov_matrix[(i, j)] == 1]

   # Mandatory cut! Otherwise, the model can find an infeasible solution
   function DegreeCut()
     x_value=Dict()
      for e in E
         value::Float64 = get_value(mctp.optimizer, x[e])
         x_value[e]=value
         if x_value[e]<=0.99 && x_value[e]>=0.01
            return
         end
      end
      for o in O ∪ M
         aux=0.0
         e=(0,0)
         for (i,j) in δ′(data, Set([o]))
            aux+=x_value[(i,j)]
            if x_value[(i,j)]>=0.999
               e=(i,j)
            end
         end
         if aux>=0.99 && aux<=1.001
            A=[]
            B=[]
            for a in δ′(data, Set([o]))
               if a!=e
                  push!(A, x[a])
                  push!(B, -1.0)
               end
            end
            push!(A, x[e])
            push!(B, 1.0)
            add_dynamic_constr!(mctp.optimizer, A, B, <=, 0.0, "degreecut")
         end 
      end
   end

   # Cover Capacity Cut (CCC)
   function CoverCapCutCallback()
      x_value = Dict()
      for e in E
         value::Float64 = get_value(mctp.optimizer, x[e])
         if value > 0.0001
            x_value[e]=value
         end
      end
      
      F = vcat([0], O, M)
      Ē = [e for e in keys(x_value)]
      data.Time+=@elapsed begin
         for k in 0:data.Max_k
            model = Model(solver = CplexSolver(CPX_PARAM_SCRIND=0))
            @variables(model, 
            begin
               x_[e in Ē], Bin
               z[i in F], Bin
               y[i in C]>=0
            end)
      
            @objective(model, Min, sum(x_value[e] * x_[e] for e in Ē))
            @constraint(model, sum(y[i] for i in C) + sum(z[j] for j in M) >= k*data.p+0.01)
            @constraint(model, [e in Ē], x_[e]>=z[e[1]] - z[e[2]])
            @constraint(model, [e in Ē], x_[e]>=z[e[2]] - z[e[1]])
            @constraint(model, [i in O], sum(y[j] for j in C if (i,j) in cov) <= z[i])
            @constraint(model, z[0] == 0)

            solve(model)
            fo = getobjectivevalue(model)
            zs = [i for i in vcat(O, M) if getvalue(z[i])>=0.5]
            sum = 0.0
            for (i,j) in Ē
               if (i ∈ zs && j ∉ zs) || (i ∉ zs && j ∈ zs)
                  sum+=x_value[(i,j)]
               end
            end
            
            if fo < 2*(k+1) - 0.001
               A=[]
               B=[]
               for (i,j) in E
                  if (i in zs && !(j in zs)) || (!(i in zs) && j in zs)
                     push!(A, x[(i,j)])
                     push!(B, 1.0)
                  end
               end
               add_dynamic_constr!(mctp.optimizer, A, B, >=, 2*(k+1), "CCC")
               data.Num_cut+=1
            end
         end
      end
   end

   # Add the cut callback to the model
   add_cut_callback!(mctp, DegreeCut, "degreecut")
   add_cut_callback!(mctp, CoverCapCutCallback, "CCC")
   
   return (mctp, x)
end
