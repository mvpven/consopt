#= # # # # # # # # # # # # # # # # # # # #
# 
# Constellation Design Optimization
#  using MGEO
# 
# # # # # # # # # # # # # # # # # # # # =#

#=== 

Case real: 2 estações e 5 PCD nos extremos, 14 satélites

===#



#incluir módulos e arquivos externos
using SatelliteToolbox
using Dates
using PyPlot
using Statistics
using MGEO
using JLD2
using LinearAlgebra

include("functions_opt.jl")
include("plot_functions.jl")


### start Data
### Ground stations and constellation design restrictions are defined here
### Argument to start_data is the initial number os satellites in the constellation.
t, epoca, pcd_list, rec_list, var_min, var_max, var_names, var_bits, nsat, pECI, pECEF, passo = start_data(14);

### calculations

function f_obj(vars) # min tgap, min a (semi-axis major)
    constellation = gen_sats_cons(vars, epoca, nsat);
    tgap, valid = satcon_optm(constellation, t, pcd_list, rec_list, pECI, pECEF; tstep = passo);
    f1 = nothing;
    f2 = nothing;
    if valid
        f1 = tgap;
        f2 = mean_a(constellation);
    end
    f = [f1,f2]

    return (valid, f)
end

function compute_pf()
    # Configuração das variáveis de projeto.
    dv = conf_design_vars(MGEO_Var(),           # ..................... Vamos usar o MGEOvar
                                   var_bits,    # .. Bits de resolução de variáveis de projeto
                                   var_min,     # ................... Valor mínimo da variável
                                   var_max,     # ................... Valor máximo da variável
                                   var_names)   # ......................... Nome das variáveis

    # Configuração do MGEO.
    mgeod = conf_mgeo(2,        # .................. Número de funções objetivo
                      1.5,      # .................. Parâmetro τ 1.5
                      5*400,    # .................. Número máximo de iterações
                      5,        # .................. Número de rodadas independentes
                      dv)       # .................. Variáveis de projeto

    # Rodar o MGEO para obter a fronteira de Pareto.
    pf = mgeo_run(mgeod, # ................................ Configuração do MGEO
                  f_obj, # .................................... Funções objetivo
                  true)  # ....................... Imprimir informações de debug

    sort_pareto!(pf,1)

    return pf
end

# Plots Pareto Frontier using  PyPlot.jl.
function plot_pf(pf)
    # Obter os valores de f1 e f2 em um grande vetor para ser plotado.
    f1 = map(x->x.f[1], pf)
    f2 = map(x->x.f[2], pf)

    # Plotar.
    plot(f1,f2, linestyle = "none", marker = "o")
    title("Pareto Frontier")
    xlabel("Fobj " * L"f_1" * " tgap [min]")
    ylabel("Fobj " * L"f_2" * " a [m]")
    grid()

end

# Execute run_opt() to run optimization. 
# testname can be changed to be used as argument to run run_opt() for several cases
function run_opt()
    testname = "application_case_sat_14"
    t_start = replace(string(DateTime(now())), ":" => "-");
    pf = compute_pf();
    t_end = replace(string(DateTime(now())), ":" => "-");
    #plot_pf(pf)
    #gcf()
    @save testname*".jdl2" pf
    println("")
    println("Começou em:    $t_start")
    println("Finalizou em:  $t_end")
    #send_email(t_start, t_end);
end

function print_data(pf)
    for point in pf
        println(point.f[1],"        ",point.f[2])
    end
end

function print_vars(pf)
    for point in pf
        println(point.vars)
    end
end