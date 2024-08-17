using JuMP
import Unicode

mutable struct Vertex
   id_vertex::Int
   pos_x::Float64
   pos_y::Float64
end

# Undirected graph
mutable struct InputGraph
   V′::Array{Vertex} # set of vertices (access with id_vertex + 1)
   E1::Array{Tuple{Int64,Int64}} # Complete set of edges
   E2::Array{Tuple{Int64,Int64}} # Set of edges not considering customers
   cost::Dict{Tuple{Int64,Int64},Float64} # cost for each edge
   cov_matrix::Dict{Tuple{Int64,Int64},Int64} # cov_matrix(ij) = 1 if the optinal node i covers the customer node j 
end

mutable struct DataCVRP
   G′::InputGraph
   M::Array{Int64} # set of mandatory vertices
   O::Array{Int64} # set of optinal vertices
   C::Array{Int64} # set of customeres
   p::Int64 #number of vehicles
   radius::Vector{Float64}
   cov::Dict{Int64,Int64} #number of times that a given customer has to be covered
   L::Float64 #Maximum time of a tour
   Time::Float64 # tempo separação CCC
   Num_cut::Int64 # quantidade cortes CCC
   Time_TwoPath::Float64 #tempo separação corte TripleCustomerQ
   NumCuts_TwoPath::Int64 # quantidade de cortes TripleCustomerQ
   Max_k::Int64 # maximum k for cut
   Name::String
end

# Read the instance informations

function read_data(app::Dict{String,Any})
   str = Unicode.normalize(read(app["instance"], String); stripcc=true)
   breaks_in = [' '; ':'; '\n']
   aux = split(str, breaks_in; limit=0, keepempty=false)
   instance_name = split(basename(app["instance"]), ".")[1]

   G′ = InputGraph([],[],[],Dict(),Dict())
   data = DataCVRP(G′, 
   [], [], [], 0, [], 
   Dict(), 10000000.0, 
   0.0, 0, 0.0, 0, 1, 
   instance_name)

   O = parse(Int, aux[2])
   C = parse(Int, aux[3])
   M = parse(Int, aux[4])-1
   data.p = parse(Int, aux[5])

   data.M = [i for i=1:M]
   data.O = [length(data.M) + i for i=1:O]
   data.C = [length(data.M) + length(data.O) + i for i=1:C]
   pos = 9
   for i=0:M+O-1
      for j=i+1:M+O
         e = (i,j)
         push!(G′.E2, e)
         data.G′.cost[e] = parse(Float64, aux[pos])
         pos += 3
      end
   end

   pos = pos - 1

   for i in data.O
      for j in data.C
         G′.cov_matrix[(i,j)] = parse(Int, aux[pos])
         
         pos = pos + 1
      end
      pos = pos + 1
   end
   pos = pos - 1


   for c in data.C
      if app["model"] == 2
         data.cov[c] = parse(Int, aux[pos])
         pos+=1
      else
         data.cov[c] = 1
      end
   end
   
   for i in 1:length(aux)
      if contains(aux[i], "NODE_COORD_SECTION")
         j = i+1
         while aux[j] != "EOF" 
            v = Vertex(0, 0, 0)
            v.id_vertex = parse(Int, aux[j])-1 # depot is forced to be 0, fisrt customer to be 1, and so on
            v.pos_x = parse(Float64, aux[j+1])
            v.pos_y = parse(Float64, aux[j+2])
            push!(G′.V′, v) # add v in the vertex array
            j+=3
         end
      end
   end


   if instance_name[1] == 't'
      data.radius = [i for i in data.O]
   else
      if instance_name[1] != 'X'
         radius = c_calculus(data)
         data.radius = [radius for i=1:length(data.O)]
      else
         radius = 1000/sqrt(O+M+1)
         data.radius = [radius for i=1:length(data.O)]

      end
   end

   cov = [(i, j) for (i, j) in keys(data.G′.cov_matrix) if data.G′.cov_matrix[(i, j)] == 1]

   # MIP to find the minimum number of vehicles to solve the problem

   function find_k(cov::Vector{Tuple{Int, Int}}, C::Vector{Int}, M::Vector{Int}, O::Vector{Int})
      model = Model(solver = CplexSolver(CPX_PARAM_SCRIND=0))
      @variables(model, begin
         x[i in C]>=0
         k, Int
      end
      )
      @objective(model, Max, k)
      @constraint(model, sum(x[i] for i in C) >= k*data.p + 0.01 - length(M))
      @constraint(model, [i in O], sum(x[j] for j in C if (i,j) in cov) <= 1)
      solve(model)
      return getobjectivevalue(model)
   end

   data.Max_k = round(Int, find_k(cov, data.C, data.M, data.O))


   max_dist = 0.0
   V = [i-1 for i=1:1 + length(data.O) + length(data.M)]
   for i in V
      if i != 0
         if c(data,(0, i)) > max_dist
            max_dist = c(data,(0, i))
         end
      end
   end
   β = 2*max_dist
   data.L = app["maxL"]+β

   return data
end


# Auxiliary functions


# Euclidian distance
function distance(data::DataCVRP, arc::Tuple{Int64, Int64})
   u, v = arc
   vertices = data.G′.V′ 
   
   x_sq = (vertices[v+1].pos_x - vertices[u+1].pos_x)^2
   y_sq = (vertices[v+1].pos_y - vertices[u+1].pos_y)^2
   #return floor(sqrt(x_sq + y_sq)+0.5)
   return sqrt(x_sq + y_sq)
end


function c_calculus(data::DataCVRP)
   aux1 = 0.0
   for v in data.O
      if !(v in data.M)
         aux = 10000000000.0
         for c in data.C
            e = (v,c)
            if distance(data,e)>0.0001 && aux > distance(data,e)
               aux = distance(data,e)
            end
         end
         if aux1 < aux
            aux1 = aux
         end
      end
   end



   aux2 = 0.0
   for c in data.C
      aux = 10000000.0
      first = 0
      for v in data.O
         if !(v in data.M)
            e = (v,c)
            if distance(data,e)>0.0001 && distance(data,e) < aux
               first = v
               aux = distance(data,e)
            end
         end
      end

      aux = 10000000.0
      second = 0
      for v in data.O
         if !(v in data.M) && v != first
            e = (v,c)
            if distance(data,e)>0.0001 && distance(data,e) < aux
               second = v
               aux = distance(data,e)
            end
         end
      end

      if aux2 < aux
         aux2 = aux
      end
   end

   if aux2 > aux1
      return aux2
   else
      return aux1
   end
end



contains(p, s) = findnext(s, p, 1) != nothing

function c(data,a) 
   if a[1] == a[2]
      return 0.0
   elseif a[1] < a[2]
      return data.G′.cost[a]
   else
      return data.G′.cost[(a[2],a[1])]
   end
end

# return incident edges of i not involving customers
function δ′(data::DataCVRP, S::Set{Int})
   incident_edges = Tuple{Int, Int}[]
   
   for (i,j) in data.G′.E2
      if i ∈ S || j ∈ S
         push!(incident_edges,(i,j))
      end
   end
   return incident_edges
end

function δ(data::DataCVRP, S)
   incident_edges = Tuple{Int, Int}[]
   for (i,j) in data.G′.E2
      if (i ∈ S && j ∉ S) || (i ∉ S && j ∈ S )
         push!(incident_edges,(i,j))
      end
   end
   return incident_edges
end

#Set of optional points that cover a given customer
function cover(data::DataCVRP, c::Integer)
   cobre = Int[]
   for i in data.O
      if data.G′.cov_matrix[(i,c)] == 1
         push!(cobre,i)
      end
   end
   return Set(cobre)
end

function ϕ(data::DataCVRP, c::Integer)
   cobre = Int[]
   for i in data.O
      if data.G′.cov_matrix[(i,c)] == 1
         push!(cobre,i)
      end
   end
   return BitSet(cobre)
end

# Set of edges e such that at least one end point of e cover the given customer c 
function cover_edges(data::DataCVRP, c::Integer)
   cobre = Tuple{Int,Int}[]
   for e in data.G′.E2
      if !(e in cobre) && ((e[1] in cover(data,c))||(e[2] in cover(data,c)))
         push!(cobre, e)
      end
   end
   return cobre
end

function Δ(data::DataCVRP, e)
   d=Int[]
   for c in data.C
      if (e[1] in data.O && data.G′.cov_matrix[(e[1],c)]==1) || (e[2] in data.O && data.G′.cov_matrix[(e[2],c)]==1)
         push!(d, c)
      end
   end
   return d
end

function σ(data::DataCVRP, i, e)
   sigma=1
   if e[1] ∈ data.O && e[2] ∈ data.O
      if e[1] ∈ cover(data, i) && e[2] ∈ cover(data, i)
         if data.cov[i] <= 1
            sigma = 0
         else
            sigma = 2
         end
      end
   end
   return sigma
end

function Σ(data, i)
   Sigma=[]
   for e in data.G′.E2
      if e[1] ∈ data.O && e[2] ∈ data.O
         if (e[1] ∉ cover(data, i) && e[2] ∈ cover(data, i)) || (e[1] ∈ cover(data, i) && e[2] ∉ cover(data, i))
            push!(Sigma, e)
         end
      end
   end
   return Sigma
end

function cover_by_optional(data::DataCVRP, o::Integer)
   cust = Int[]
   for c in data.C
      if data.G′.cov_matrix[(o,c)] == 1
         push!(cust,c)
      end
   end
   return cust
end