if !("." in LOAD_PATH)
    push!(LOAD_PATH, ".")
end
import Coverage
import Plots

# Packages to monitor
import AbstractFFTs
import AbstractTrees
import Adapt
import ArgTools
import ArnoldiMethod
import ArrayInterface
import ArrayLayouts
import Artifacts
import BFloat16s
import BandedMatrices
import Bijections
import BinaryProvider
import BitTwiddlingConvenienceFunctions
import BlockArrays
import BoundaryValueDiffEq
import CEnum
import CFTime
import CLIMAParameters
import CPUSummary
import CSTParser
import CUDA
import Cassette
import ChainRulesCore
import CloseOpenIntervals
import CodeTracking
import Combinatorics
import CommonMark
import CommonSolve
import CommonSubexpressions
import Compat
import CompositeTypes
import ConstructionBase
import CpuId
import Crayons
import DEDataArrays
import DataAPI
import DataStructures
import DataValueInterfaces
import DelayDiffEq
import DiffEqBase
import DiffEqCallbacks
import DiffEqDiffTools
import DiffEqFinancial
import DiffEqJump
import DiffEqNoiseProcess
import DiffEqPhysics
import DiffResults
import DiffRules
import Distances
import Distributions
import DomainSets
import Downloads
import DynamicPolynomials
import EllipsisNotation
import ExponentialUtilities
import ExprTools
import FastBroadcast
import FastClosures
import FillArrays
import FiniteDiff
import Formatting
import ForwardDiff
import FunctionWrappers
import GPUArrays
import GPUCompiler
import GaussQuadrature
import Highlights
import Hwloc
import IfElse
import Inflate
import IntervalSets
import IterativeSolvers
import IteratorInterfaceExtensions
import JLLWrappers
import JSON
import JuliaFormatter
import JuliaInterpreter
import KernelAbstractions
import LLVM
import LaTeXStrings
import LabelledArrays
import Latexify
import LayoutPointers
import LazyArtifacts
import LeftChildRightSiblingTrees
import LibCURL
import LightGraphs
import LineSearches
import LogExpFunctions
import LoopVectorization
import MacroTools
import ManualMemory
import Memoize
import Missings
import ModelingToolkit
import MuladdMacro
import MultiScaleArrays
import MultivariatePolynomials
import MutableArithmetics
import NLSolversBase
import NLsolve
import NaNMath
import NetworkOptions
import NonlinearSolve
import NonlinearSolvers
import OffsetArrays
import OrderedCollections
import OrdinaryDiffEq
import PDMats
import ParameterizedFunctions
import Parameters
import Parsers
import PoissonRandom
import Polyester
import PolyesterWeave
import Polynomials
import PreallocationTools
import ProgressLogging
import QuadGK
import RandomNumbers
import RecipesBase
import RecursiveArrayTools
import RecursiveFactorization
import Reexport
import Requires
import ResettableStacks
import Rmath
import RootSolvers
import Roots
import RuntimeGeneratedFunctions
import SIMDTypes
import SLEEFPirates
import SafeTestsets
import SciMLBase
import Scratch
import Setfield
import SimpleTraits
import SortingAlgorithms
import SparseDiffTools
import SpecialFunctions
import Static
import StaticArrays
import StatsBase
import StatsFuns
import SteadyStateDiffEq
import StochasticDiffEq
import StrideArraysCore
import Sundials
import SymbolicUtils
import Symbolics
import TableTraits
import Tables
import TermInterface
import TextWrap
import Thermodynamics
import ThreadingUtilities
import TimerOutputs
import Tokenize
import TreeViews
import URIParser
import UnPack
import Unitful
import VectorizationBase
import VertexSafeGraphs
import ZygoteRules
import ClimaCore
import CloudMicrophysics
import Debugger
import Dierckx
import DifferentialEquations
import FastGaussQuadrature
import Glob
import LambertW
import NCDatasets
import PrettyTables
import SurfaceFluxes

mod_dir(x) = dirname(dirname(pathof(x)))
all_dirs_to_monitor = [
    ".",
    "/central/software/julia/1.7.0/share/julia/base/",
    mod_dir(AbstractFFTs),
    mod_dir(AbstractTrees),
    mod_dir(Adapt),
    mod_dir(ArgTools),
    mod_dir(ArnoldiMethod),
    mod_dir(ArrayInterface),
    mod_dir(ArrayLayouts),
    mod_dir(Artifacts),
    mod_dir(BFloat16s),
    mod_dir(BandedMatrices),
    mod_dir(Bijections),
    mod_dir(BinaryProvider),
    mod_dir(BitTwiddlingConvenienceFunctions),
    mod_dir(BlockArrays),
    mod_dir(BoundaryValueDiffEq),
    mod_dir(CEnum),
    mod_dir(CFTime),
    mod_dir(CLIMAParameters),
    mod_dir(CPUSummary),
    mod_dir(CSTParser),
    mod_dir(CUDA),
    mod_dir(Cassette),
    mod_dir(ChainRulesCore),
    mod_dir(CloseOpenIntervals),
    mod_dir(CodeTracking),
    mod_dir(Combinatorics),
    mod_dir(CommonMark),
    mod_dir(CommonSolve),
    mod_dir(CommonSubexpressions),
    mod_dir(Compat),
    mod_dir(CompositeTypes),
    mod_dir(ConstructionBase),
    mod_dir(CpuId),
    mod_dir(Crayons),
    mod_dir(DEDataArrays),
    mod_dir(DataAPI),
    mod_dir(DataStructures),
    mod_dir(DataValueInterfaces),
    mod_dir(DelayDiffEq),
    mod_dir(DiffEqBase),
    mod_dir(DiffEqCallbacks),
    mod_dir(DiffEqDiffTools),
    mod_dir(DiffEqFinancial),
    mod_dir(DiffEqJump),
    mod_dir(DiffEqNoiseProcess),
    mod_dir(DiffEqPhysics),
    mod_dir(DiffResults),
    mod_dir(DiffRules),
    mod_dir(Distances),
    mod_dir(Distributions),
    mod_dir(DomainSets),
    mod_dir(Downloads),
    mod_dir(DynamicPolynomials),
    mod_dir(EllipsisNotation),
    mod_dir(ExponentialUtilities),
    mod_dir(ExprTools),
    mod_dir(FastBroadcast),
    mod_dir(FastClosures),
    mod_dir(FillArrays),
    mod_dir(FiniteDiff),
    mod_dir(Formatting),
    mod_dir(ForwardDiff),
    mod_dir(FunctionWrappers),
    mod_dir(GPUArrays),
    mod_dir(GPUCompiler),
    mod_dir(GaussQuadrature),
    mod_dir(Highlights),
    mod_dir(Hwloc),
    mod_dir(IfElse),
    mod_dir(Inflate),
    mod_dir(IntervalSets),
    mod_dir(IterativeSolvers),
    mod_dir(IteratorInterfaceExtensions),
    mod_dir(JLLWrappers),
    mod_dir(JSON),
    mod_dir(JuliaFormatter),
    mod_dir(JuliaInterpreter),
    mod_dir(KernelAbstractions),
    mod_dir(LLVM),
    mod_dir(LaTeXStrings),
    mod_dir(LabelledArrays),
    mod_dir(Latexify),
    mod_dir(LayoutPointers),
    mod_dir(LazyArtifacts),
    mod_dir(LeftChildRightSiblingTrees),
    mod_dir(LibCURL),
    mod_dir(LightGraphs),
    mod_dir(LineSearches),
    mod_dir(LogExpFunctions),
    mod_dir(LoopVectorization),
    mod_dir(MacroTools),
    mod_dir(ManualMemory),
    mod_dir(Memoize),
    mod_dir(Missings),
    mod_dir(ModelingToolkit),
    mod_dir(MuladdMacro),
    mod_dir(MultiScaleArrays),
    mod_dir(MultivariatePolynomials),
    mod_dir(MutableArithmetics),
    mod_dir(NLSolversBase),
    mod_dir(NLsolve),
    mod_dir(NaNMath),
    mod_dir(NetworkOptions),
    mod_dir(NonlinearSolve),
    mod_dir(NonlinearSolvers),
    mod_dir(OffsetArrays),
    mod_dir(OrderedCollections),
    mod_dir(OrdinaryDiffEq),
    mod_dir(PDMats),
    mod_dir(ParameterizedFunctions),
    mod_dir(Parameters),
    mod_dir(Parsers),
    mod_dir(PoissonRandom),
    mod_dir(Polyester),
    mod_dir(PolyesterWeave),
    mod_dir(Polynomials),
    mod_dir(PreallocationTools),
    mod_dir(ProgressLogging),
    mod_dir(QuadGK),
    mod_dir(RandomNumbers),
    mod_dir(RecipesBase),
    mod_dir(RecursiveArrayTools),
    mod_dir(RecursiveFactorization),
    mod_dir(Reexport),
    mod_dir(Requires),
    mod_dir(ResettableStacks),
    mod_dir(Rmath),
    mod_dir(RootSolvers),
    mod_dir(Roots),
    mod_dir(RuntimeGeneratedFunctions),
    mod_dir(SIMDTypes),
    mod_dir(SLEEFPirates),
    mod_dir(SafeTestsets),
    mod_dir(SciMLBase),
    mod_dir(Scratch),
    mod_dir(Setfield),
    mod_dir(SimpleTraits),
    mod_dir(SortingAlgorithms),
    mod_dir(SparseDiffTools),
    mod_dir(SpecialFunctions),
    mod_dir(Static),
    mod_dir(StaticArrays),
    mod_dir(StatsBase),
    mod_dir(StatsFuns),
    mod_dir(SteadyStateDiffEq),
    mod_dir(StochasticDiffEq),
    mod_dir(StrideArraysCore),
    mod_dir(Sundials),
    mod_dir(SymbolicUtils),
    mod_dir(Symbolics),
    mod_dir(TableTraits),
    mod_dir(Tables),
    mod_dir(TermInterface),
    mod_dir(TextWrap),
    mod_dir(Thermodynamics),
    mod_dir(ThreadingUtilities),
    mod_dir(TimerOutputs),
    mod_dir(Tokenize),
    mod_dir(TreeViews),
    mod_dir(URIParser),
    mod_dir(UnPack),
    mod_dir(Unitful),
    mod_dir(VectorizationBase),
    mod_dir(VertexSafeGraphs),
    mod_dir(ZygoteRules),
    mod_dir(ClimaCore),
    mod_dir(CloudMicrophysics),
    mod_dir(Debugger),
    mod_dir(Dierckx),
    mod_dir(DifferentialEquations),
    mod_dir(FastGaussQuadrature),
    mod_dir(Glob),
    mod_dir(LambertW),
    mod_dir(NCDatasets),
    mod_dir(PrettyTables),
    mod_dir(SurfaceFluxes),
]

# https://github.com/jheinen/GR.jl/issues/278#issuecomment-587090846
ENV["GKSwstype"] = "nul"

all_cases = [
    # "ARM_SGP",
    # "Bomex",
    # "DryBubble",
    # "DYCOMS_RF01",
    # "GABLS",
    # "GATE_III",
    # "life_cycle_Tan2018",
    # "Nieuwstadt",
    "Rico",
    # "Soares",
    # "SP",
    "TRMM_LBA",
    "LES_driven_SCM",
]

filter!(x -> x ≠ "GATE_III", all_cases) # no mse tables for GATE_III
filter!(x -> x ≠ "SP", all_cases) # not currently running SP
allocs = Dict()
for case in all_cases
    ENV["ALLOCATION_CASE_NAME"] = case
    run(`julia --project=test/ --track-allocation=all perf/alloc_per_case.jl`)

    allocs[case] = Coverage.analyze_malloc(all_dirs_to_monitor)

    # Clean up files
    for d in all_dirs_to_monitor
        all_files = [joinpath(root, f) for (root, dirs, files) in Base.Filesystem.walkdir(d) for f in files]
        all_mem_files = filter(x -> endswith(x, ".mem"), all_files)
        for f in all_mem_files
            rm(f)
        end
    end
end

@info "Post-processing allocations"

function plot_allocs(case_name, allocs_per_case, n_unique_bytes)
    p = Plots.plot()
    @info "Allocations for $case_name"

    function filename_only(fn)
        fn = first(split(fn, ".jl")) * ".jl"
        splitby = "central/scratch/climaci/turbulenceconvection-ci/depot/cpu/packages/"
        if occursin(splitby, fn)
            fn = last(split(fn, splitby))
        end
        return fn
    end
    function compile_tc(fn, linenumber)
        c1 = endswith(filename_only(fn), "TurbulenceConvection.jl")
        c2 = linenumber == 1
        return c1 && c2
    end

    filter!(x -> x.bytes ≠ 0, allocs_per_case)
    filter!(x -> !compile_tc(x.filename, x.linenumber), allocs_per_case)

    for alloc in allocs_per_case
        println(alloc)
    end
    println("Number of allocating sites: $(length(allocs_per_case))")
    case_bytes = getproperty.(allocs_per_case, :bytes)[end:-1:1]
    case_filename = getproperty.(allocs_per_case, :filename)[end:-1:1]
    case_linenumber = getproperty.(allocs_per_case, :linenumber)[end:-1:1]
    all_bytes = Int[]
    filenames = String[]
    linenumbers = Int[]
    loc_ids = String[]
    for (bytes, filename, linenumber) in zip(case_bytes, case_filename, case_linenumber)
        compile_tc(filename, linenumber) && continue # Skip loading module
        loc_id = "$(filename_only(filename))" * "$linenumber"
        if !(bytes in all_bytes) && !(loc_id in loc_ids)
            push!(all_bytes, bytes)
            push!(filenames, filename)
            push!(linenumbers, linenumber)
            push!(loc_ids, loc_id)
            if length(all_bytes) ≥ n_unique_bytes
                break
            end
        end
    end

    all_bytes = all_bytes ./ 10^3
    max_bytes = maximum(all_bytes)
    @info "$case_name: $all_bytes"
    xtick_name(filename, linenumber) = "$filename, line number: $linenumber"
    markershape = (:square, :hexagon, :circle, :star, :utriangle, :dtriangle)
    for (bytes, filename, linenumber) in zip(all_bytes, filenames, linenumbers)
        Plots.plot!(
            [0],
            [bytes];
            seriestype = :scatter,
            label = xtick_name(filename_only(filename), linenumber),
            markershape = markershape[1],
            markersize = 1 + bytes / max_bytes * 10,
        )
        markershape = (markershape[end], markershape[1:(end - 1)]...)
    end
    p1 = Plots.plot!(ylabel = "Allocations (KB)", title = case_name)
    subset_allocs_per_case = allocs_per_case[end:-1:(end - 100)]
    p2 = Plots.plot(
        1:length(subset_allocs_per_case),
        getproperty.(subset_allocs_per_case, :bytes) ./ 1000;
        xlabel = "i-th allocating line (truncated and sorted)",
        ylabel = "Allocations (KB)",
        markershape = :circle,
    )
    Plots.plot(p1, p2, layout = Plots.grid(2, 1))
    Plots.savefig(joinpath(folder, "allocations_$case_name.png"))
end

folder = "perf/allocations_output"
mkpath(folder)

@info "Allocated bytes for single tendency per case:"
for case in all_cases
    plot_allocs(case, allocs[case], 10)
end
