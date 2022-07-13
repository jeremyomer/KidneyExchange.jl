function create_model(time_limit::Float64, optimizer::String, is_integer::Bool = false, verbose::Bool = false)
  if optimizer=="CPLEX"
    model = Model(CPLEX.Optimizer)
    set_optimizer_attribute(model, "CPX_PARAM_THREADS", 1)
    set_optimizer_attribute(model, "CPX_PARAM_TILIM", max(0.0,time_limit))
    if verbose
      set_optimizer_attribute(model, "CPX_PARAM_MIPDISPLAY", 2)
    else
      set_optimizer_attribute(model, "CPX_PARAM_SCRIND", 0)
    end
  elseif optimizer=="Gurobi"
    model = Model(Gurobi.Optimizer)
    set_optimizer_attribute(model, "OutputFlag", Int(verbose))
    set_optimizer_attribute(model, "TimeLimit", max(0.0,time_limit))
    set_optimizer_attribute(model, "Threads", 1)
  elseif optimizer=="GLPK"
    model = Model(GLPK.Optimizer)
    set_optimizer_attribute(model, "tm_lim", round(Int, time_limit) * 1_000)
    if verbose set_optimizer_attribute(model, "msg_lev", 3) end
  elseif optimizer=="GLPK-Cbc"
    if is_integer
      model = Model(Cbc.Optimizer)
      set_optimizer_attribute(model, "seconds", time_limit)
      set_optimizer_attribute(model, "logLevel", Int(verbose))
    else
      model = Model(GLPK.Optimizer)
      set_optimizer_attribute(model, "tm_lim", round(Int, time_limit) * 1_000)
      if verbose set_optimizer_attribute(model, "msg_lev", 3) end
    end
  elseif optimizer=="Clp"
    if is_integer
      model = Model(Cbc.Optimizer)
      set_optimizer_attribute(model, "seconds", time_limit)
      set_optimizer_attribute(model, "logLevel", Int(verbose))
    else
      model = Model(Clp.Optimizer)
      set_optimizer_attribute(model, "MaximumSeconds", time_limit)
      set_optimizer_attribute(model, "LogLevel", Int(verbose))
    end
  elseif optimizer=="Cbc"
    if is_integer
      model = Model(Cbc.Optimizer)
      set_optimizer_attribute(model, "seconds", time_limit)
      set_optimizer_attribute(model, "logLevel", Int(verbose))
    else
      model = Model(Clp.Optimizer)
      set_optimizer_attribute(model, "MaximumSeconds", time_limit)
      set_optimizer_attribute(model, "LogLevel", Int(verbose))
    end
  else
    println(optimizer)
    error("The chosen optimizer is unknown!")
  end
  return model
end

function set_time_limit(model::Model, time_limit::Float64, optimizer::String)
  if optimizer=="CPLEX"
    set_optimizer_attribute(model, "CPX_PARAM_TILIM", max(0.0,time_limit))
  elseif optimizer=="Gurobi"
    set_optimizer_attribute(model, "TimeLimit", max(0.0,time_limit))
  elseif optimizer=="GLPK"
    set_optimizer_attribute(model, "tm_lim", round(Int, time_limit) * 1_000)
  elseif optimizer=="GLPK-Cbc"
    set_optimizer_attribute(model, "seconds", time_limit)
  elseif optimizer=="Clp"
      set_optimizer_attribute(model, "seconds", time_limit)
  elseif optimizer=="Cbc"
      set_optimizer_attribute(model, "seconds", time_limit)
  else
    println(optimizer)
    error("The chosen optimizer is unknown!")
  end
end
