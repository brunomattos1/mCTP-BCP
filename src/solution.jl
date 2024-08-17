mutable struct Solution
   cost::Union{Int,Float64}
   routes::Array{Array{Int}}
end

# build Solution from the variables x
function getsolution(data::DataCVRP, optimizer::VrpOptimizer, x, objval, app::Dict{String,Any})
   E = data.G′.E2
   O = data.O
   routes = []

   cont = 1
   for path_id in 1:get_number_of_positive_paths(optimizer)
      route_edeges, route, visited, visits =[], [], Dict(), Dict()
      cont += 1
      for e in E
         visits[e] = 0
         visited[e] = 0
         if get_value(optimizer, x[e], path_id) > 0.5
            push!(route_edeges,e)
            visited[e] = floor(get_value(optimizer, x[e], path_id)+0.5)
         end
      end
      current = 0
      total_e = 0
      for e in route_edeges
         if e[1] == 0
            visits[e] = visits[e] + 1
            route = [0, e[2]]
            current = e[2]
            break
         end
      end
      for e in route_edeges
         total_e += visited[e]
      end
      number_e = 1

      while current != 0
         E1 = []
         for e in route_edeges
            if visits[e] < visited[e] && (e[1] == current || e[2] == current)
               if number_e < total_e - 1 && e[1] == 0
                  push!(E1, e)
               else
                  number_e += 1
                  visits[e] = visits[e] + 1
                  if e[1] == current
                     current = e[2]
                  else 
                     current = e[1]
                  end
                  push!(route, current)
               end
            end
         end
         for e in E1
            if visits[e] < visited[e] && (e[1] == current || e[2] == current)
               visits[e] = visits[e] + 1
               if e[1] == current
                  current = e[2]
               else 
                  current = e[1]
               end
               push!(route, current)
            end
         end
      end
      push!(routes, route)
   end
   return Solution(objval, routes)
end

function getsolution2(data, optimizer, x, objval, app)
   E = data.G′.E2
   O = data.O
   M = data.M
   routes = []
   T = [i for i = 1:length(M) + length(O) - app["veh"] + 1]
   aux = []
   for path_id in 1:get_number_of_positive_paths(optimizer)
      push!(routes, [])
      push!(aux, [])
      for e in E, t in T
         if get_value(optimizer, x[e,t], path_id) > 0.5
            push!(routes[path_id], e)
            push!(aux[path_id], e[1], e[2])
         end
      end
   end
   #@show aux[1]
   for r=1:length(aux)
      for i=2:length(aux[r])
         cont = 1
         for j=2:length(aux[r])
            if i != j
               if aux[r][i] == aux[r][j]
                  cont += 1
               end
            end
         end
         if cont == 1
            push!(routes[r], (0, aux[r][i]))
         end
      end
   end


   costs = zeros(length(routes))
   for r=1:length(routes)
      for e in routes[r]
         costs[r] += c(data, e)
      end
   end

   routes_ = []
   for r = 1:length(routes)
      push!(routes_, [])
      append!(routes_[r], routes[r][1])
      while length(routes_[r]) < length(routes[r])
         for e in routes[r]
            if 0 ∉ e
               if length(e ∩ routes_[r][end]) == 1
                  append!(routes_[r], setdiff(e, routes_[r][end]))
               end
            end
         end
      end
      push!(routes_[r], 0)
   end

   return Solution(sum(costs), routes_)
end
# Print the solution
function print_routes(solution)
   println("############################")
   for (i,r) in enumerate(solution.routes)
      print("Route #$i: ") 
      for j in r
         print("$j ")
      end
      println()
   end
   println("Cost: ", solution.cost)
   println("############################")
end

# Checks the feasiblity of a solution
function checksolution(data::DataCVRP, app::Dict{String,Any}, solution)
   E, C, O, M = data.G′.E2, data.C, data.O, data.M

   cost = 0.0

   number_cov = Dict()
   for c in C
      number_cov[c] = 0
   end

   visit_M = Dict()
   for m in M
      visit_M[m] = 0
   end

   visit_O = Dict()
   for r in solution.routes
      for o in O
         visit_O[o] = 0
      end
      len = 0.0
      nv = length(r)-2
      if nv > data.p
         println("error:")
         @show r
         @show nv
         @show data.p
      end
      for i=1:length(r)-1

         if r[i+1] in O
            visit_O[r[i+1]] += 1
            for c in cover_by_optional(data,r[i+1])
               number_cov[c] += 1
            end
         end
         if r[i+1] in M
            visit_M[r[i+1]] += 1
         end
         e = (r[i],r[i+1])
         if r[i] > r[i+1]
            e = (r[i+1],r[i])
         end
         len += c(data,e)
      end
      cost += len
      for o in O
         if visit_O[o] > 1
            println("error")
            @show r
            @show o, visit_O[o]
         end
      end

      max = 0.0
      for i in vcat(O,M)
         e = (0,i)
         if c(data,e) > max
            max = c(data,e)
         end
      end
      if app["maxL"] < 10000000.0 && len > app["maxL"] + 2*max
         println("error:")
         @show r
         @show len
         @show app["maxL"] + 2*max
      end
   end

   for c in C
      if number_cov[c] < 1
         println("error:")
         @show solution.routes
         @show cover(data,c)
         @show c
         @show number_cov[c]
         @show data.cov[c]
      end
      if number_cov[c] < data.cov[c]
         println("erro data cov")
         @show solution.routes
         @show cover(data,c)
         @show c
         @show number_cov[c]
         @show data.cov[c]
      end
   end

   for m in M
      if visit_M[m] != 1
         println("error")
         @show solution.routes
         @show m, visit_M[m]
      end
   end

   if cost < solution.cost - 0.0001 || cost > solution.cost + 0.0001
      println("error ")
      @show cost
      @show solution.cost
   end
   #=
   dim, Q = dimension(data), veh_capacity(data)
   visits = [0 for i in 2:dim]
   sum_cost = 0.0
   for (i,r) in enumerate(solution.routes)
      sum_demand, prev = 0.0, 0
      for j in r
         visits[j] += 1
         (visits[j] == 2) && error("Customer $j was visited more than once")
         sum_cost += distance(data, (prev,j))
         sum_demand += d(data, j)
         prev = j
      end
      sum_cost += distance(data, (prev,0))
      (sum_demand > Q) && error("Route #$i is violating the capacity constraint. Sum of the demands is $(sum_demand) and Q is $Q")
   
   end
   !isempty(filter(a->a==0,visits)) && error("The following customers were not visited: $(filter(a->a==0,visits))")
   (abs(solution.cost-sum_cost) > 0.001) && error("Cost calculated from the routes ($sum_cost) is different from that passed as"*
   =#                                                                                              # " argument ($(solution.cost)).") 
end

# Write solution in a file
function writesolution(solpath, solution)
   open(solpath, "w") do f
      for (i,r) in enumerate(solution.routes)
         write(f, "Route #$i: ")
         for j in r
            write(f, "$j ") 
         end
         write(f, "\n")
      end
      write(f, "Cost $(solution.cost)\n")
   end
end

# Draw the solution in a .pdf
function drawsolution(tikzpath, data, solution)
   @show data.Name[1]
   V = data.G′.V′
   E, C, O, M = data.G′.E2, data.C, data.O, data.M
   #V = [i-1 for i=1:1 + length(O) + length(M)]

   desenha = Dict()
   for i in O ∪ M
      if data.Name[1] == 'X'
         desenha[i]=false
      else
         desenha[i]=true
      end
   end
   if data.Name[1] == 'X'
         
      for r in solution.routes
         for i in r
            if i in O
               desenha[i]=true
            end
         end
      end
   end


   
   open(tikzpath, "w") do f
      write(f,"\\documentclass[crop,tikz]{standalone}\n\\begin{document}\n")
      # get limits to draw
      pos_x_vals = [i.pos_x for i in V]
      pos_y_vals = [i.pos_y for i in V]
      scale_fac = 1/(max(maximum(pos_x_vals),maximum(pos_y_vals))/10)
      write(f,"\\begin{tikzpicture}[thick, scale=1, every node/.style={scale=0.3}]\n")

      for i in V
         x_plot = scale_fac*i.pos_x
         y_plot = scale_fac*i.pos_y
         if i.id_vertex == 0 # plot depot
            write(f, "\t\\node[draw, line width=0.1mm, rectangle, fill=yellow, inner sep=0.05cm, scale=1.4] (v$(i.id_vertex)) at ($(x_plot),$(y_plot)) {\\footnotesize $(i.id_vertex)};\n")
            # Uncomment to plot without vertex id
            #write(f, "\t\\node[draw, rectangle, fill=yellow, scale=1.4] (v$(i.id_vertex)) at ($(x_plot),$(y_plot)) {};\n")
         elseif i.id_vertex in vcat(M,O)
            if i.id_vertex in O
               if desenha[i.id_vertex]
                  write(f, "\t\\node[draw, line width=0.1mm, circle, fill=red, inner sep=0.05cm] (v$(i.id_vertex)) at ($(x_plot),$(y_plot)) {\\footnotesize $(i.id_vertex)};\n")
               end
            else
               write(f, "\t\\node[draw, line width=0.1mm, rectangle, fill=green, inner sep=0.05cm] (v$(i.id_vertex)) at ($(x_plot),$(y_plot)) {\\footnotesize $(i.id_vertex)};\n")
            end
               # Uncomment to plot without vertex id
            #write(f, "\t\\node[draw, circle, fill=white] (v$(i.id_vertex)) at ($(x_plot),$(y_plot)) {};\n")
         else
            if data.Name[1] != 'X' || desenha[i.id_vertex - length(O) - length(M)]==false
               write(f, "\t\\node[draw, line width=0.1mm, circle, fill=cyan, inner sep=0.05cm] (v$(i.id_vertex)) at ($(x_plot),$(y_plot)) {\\footnotesize $(i.id_vertex)};\n")
            end
         end
      end
      for r in solution.routes
         #=prev = r[1] # Uncomment (and comment below) to hide edges with the depot
         for i in r[2:end]
            e = (prev,i)
            write(f, "\t\\draw[-,line width=0.8pt] (v$(e[1])) -- (v$(e[2]));\n")
            prev = i
         end=#
         prev = 0
         for i in r
            e = (prev,i)
            if i in O
               if desenha[i]
                  radius = scale_fac*data.radius[i - length(data.M)]
                  write(f,"\t\\draw[red,dashed,-,line width=0.2pt,opacity=.4] (v$(i)) circle ($(radius));\n")
               end
            end
            #edge_style = (e[1] == 0 || e[2] == 0) ? "dashed,-,line width=0.2pt,opacity=.2" : "-,line width=0.8pt"
            edge_style =  "-,line width=0.8pt"
            write(f, "\t\\draw[$(edge_style)] (v$(e[1])) -- (v$(e[2]));\n")
            prev = i
         end
         #write(f, "\t\\draw[dashed,-,line width=0.2pt,opacity=.2] (v0) -- (v$(prev));\n") 
      end
      write(f, "\\end{tikzpicture}\n")
      write(f, "\\end{document}\n")
   end   
end

# Draw a fractional solution
function drawfrac(path, data, x_value)
   appfolder = dirname(@__FILE__)
   V = data.G′.V′
   E, C, O, M = data.G′.E2, data.C, data.O, data.M
   open(path, "w") do f
      write(f,"\\documentclass[crop,tikz]{standalone}\n\\begin{document}\n")
      # get limits to draw
      pos_x_vals = [i.pos_x for i in V]
      pos_y_vals = [i.pos_y for i in V]
      scale_fac = 1/(max(maximum(pos_x_vals),maximum(pos_y_vals))/10)
      write(f,"\\begin{tikzpicture}[thick, scale=1, every node/.style={scale=0.3}]\n")
      for i in V
         x_plot = scale_fac*i.pos_x
         y_plot = scale_fac*i.pos_y
         if i.id_vertex == 0 # plot depot
            write(f, "\t\\node[draw, line width=0.1mm, rectangle, fill=yellow, inner sep=0.05cm, scale=1.4] (v$(i.id_vertex)) at ($(x_plot),$(y_plot)) {\\footnotesize $(i.id_vertex)};\n")
            # Uncomment to plot without vertex id
            #write(f, "\t\\node[draw, rectangle, fill=yellow, scale=1.4] (v$(i.id_vertex)) at ($(x_plot),$(y_plot)) {};\n")
         elseif i.id_vertex in vcat(M,O)
            if i.id_vertex in O
               radius = scale_fac*data.radius[1]
               write(f, "\t\\node[draw, line width=0.1mm, circle, fill=red, inner sep=0.05cm] (v$(i.id_vertex)) at ($(x_plot),$(y_plot)) {\\footnotesize $(i.id_vertex)};\n")
               for e in keys(x_value)
                  if e[1] == i.id_vertex || e[2] == i.id_vertex
                     write(f,"\t\\draw[red,dashed,-,line width=0.2pt,opacity=.4] (v$(i.id_vertex)) circle ($(radius));\n")
                     break
                  end
               end
            else
               write(f, "\t\\node[draw, line width=0.1mm, rectangle, fill=green, inner sep=0.05cm] (v$(i.id_vertex)) at ($(x_plot),$(y_plot)) {\\footnotesize $(i.id_vertex)};\n")
            end
            # Uncomment to plot without vertex id
            #write(f, "\t\\node[draw, circle, fill=white] (v$(i.id_vertex)) at ($(x_plot),$(y_plot)) {};\n")
         else
            write(f, "\t\\node[draw, line width=0.1mm, circle, fill=cyan, inner sep=0.05cm] (v$(i.id_vertex)) at ($(x_plot),$(y_plot)) {\\footnotesize $(i.id_vertex)};\n")
         end
      end
      for e in keys(x_value)
         edge_style =  string("-,line width=",x_value[e],"pt")
         write(f, "\t\\draw[$(edge_style)] (v$(e[1])) -- (v$(e[2]));\n")
      end
      write(f, "\\end{tikzpicture}\n")
      write(f, "\\end{document}\n")
   end   
end

