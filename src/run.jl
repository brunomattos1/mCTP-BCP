using VrpSolver, JuMP, ArgParse, CPLEX
#using CPLEX
include("data.jl")
include("model.jl")
include("solution.jl")
#include("SparseMaxFlowMinCut.jl")

function parse_commandline(args_array::Array{String,1}, appfolder::String)
   s = ArgParseSettings(usage="##### VRPSolver #####\n\n"*
	   "  On interactive mode, call main([\"arg1\", ..., \"argn\"])", exit_after_help=false)
   @add_arg_table! s begin
      "instance"
         help = "Instance file path"
      "--cfg", "-c"
         help = "Configuration file path"
         default = "$appfolder/../config/CVRP.cfg"
      "--ub","-u"
         help = "Upper bound (primal bound)"
         arg_type = Float64
         default = 10000000.0
      "--sol","-s"
         help = "Solution file path (CVRPLIB format. "*
                  "e.g. http://vrp.atd-lab.inf.puc-rio.br/media/com_vrp/instances/E/E-n13-k4.sol)"
      "--out","-o"
         help = "Path to write the solution found"
      "--tikz","-t"
         help = "Path to write the TikZ figure of the solution found."
      "--fractex","-f"
         help = "Path to write the TikZ figure of the fractionary solution"
      "--nosolve","-n"
         help = "Does not call the VRPSolver. Only to check or draw a given solution."
         action = :store_true
      "--maxc","-p"
         help = "Maximum number of customers in a route"
         arg_type = Int
         default = 4  
      "--mand","-M"
         help = "Number of mandatory points"
         arg_type = Int
         default = 0 
      "--opt","-O"
         help = "Number of optional points"
         arg_type = Int
         default = 0 
      "--cust","-C"
         help = "Number of customers points"
         arg_type = Int
         default = 0 
      "--model","-m"
         help = "type of model: 1 - mctp; 2 - mmctp binary; 3 - mmctp without overnight; 4 - mmctp with overnight; 5 - cumulative"
         arg_type = Int
         default = 1 
      "--maxL","-q"
         help = "Maximum length for a route"
         arg_type = Float64
         default = 10000000.0
      "--veh", "-k"
         help = "max number of vehicles"
         arg_type = Int
         default = 2
      
   end
   return parse_args(args_array, s)
end

function run_cvrp(app::Dict{String,Any})
   println("Application parameters:")
   for (arg,val) in app
      println("  $arg  =>  $(repr(val))")
   end
   flush(stdout)

   instance_name = split(basename(app["instance"]), ".")[1]

   data = read_data(app)

   if app["sol"] != nothing
      sol = readsolution(app)
      app["ub"] = (sol.cost < app["ub"]) ? sol.cost : app["ub"] # update the upper bound if necessary
   end

   solution_found = false
   if !app["nosolve"]
      (model, x) = build_model(data, app)

      
      optimizer = VrpOptimizer(model, app["cfg"], instance_name)
      set_cutoff!(optimizer, app["ub"] + 0.01)

      (status, solution_found) = optimize!(optimizer)
      if solution_found
         sol = getsolution(data, optimizer, x, get_objective_value(optimizer), app)
         print_routes(sol)
      end
   end
   if solution_found || app["sol"] != nothing # Is there a solution?
      checksolution(data, app, sol)
      
      println("Info_cuts_cols: Instance & #CCC & TimeCCC & #TwoPath & TimeTwoPath")
      println("Info_cuts: $(data.Name) & $(data.Num_cut) & $(round(data.Time, digits = 2))")
      println("############################")
      
      
      if app["out"] != nothing
         writesolution(app["out"], sol)
      end
      if app["tikz"] != nothing
         drawsolution(app["tikz"], data, sol) # write tikz figure
      end
   elseif !app["nosolve"]
      if status == :Optimal
         println("Problem infeasible")
      else
         println("Solution not found")
      end
      
      println("Info_cuts_cols: Instance & #CCC & TimeCCC")
      println("Info_cuts: $(data.Name) & $(data.Num_cut) & $(round(data.Time, digits = 2))")
      println("############################")
   end
   #println("########################################################")
end

function main(args)
   appfolder = dirname(@__FILE__)
   app = parse_commandline(args, appfolder)
   isnothing(app) && return
   run_cvrp(app)
end

if isempty(ARGS)
   main(["--help"])
else
   main(ARGS)
end
