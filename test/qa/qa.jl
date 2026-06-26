using SciMLTesting, IRKGaussLegendre, JET, Test

run_qa(
    IRKGaussLegendre;
    explicit_imports = true,
    jet_kwargs = (; target_defined_modules = true),
    ei_kwargs = (;
        all_qualified_accesses_are_public = (;
            ignore = (
                :DEStats, :__solve, :allows_arbitrary_number_types, :allowscomplex,
                :forwarddiffs_model, :forwarddiffs_model_time, :isadaptive,
                :isautodifferentiable, :unwrapped_f,         # SciMLBase internals
                :Forward, :Reverse, :TwicePrecision,         # Base internals
            ),
        ),
        # `_process_verbose_param` is a non-public DiffEqBase helper.
        all_explicit_imports_are_public = (; ignore = (:_process_verbose_param,)),
    ),
    # no_implicit_imports cannot be cleanly fixed: the module does
    # `@reexport using SciMLBase`, so the SciMLBase API names are intentionally
    # implicit. Tracked in https://github.com/SciML/IRKGaussLegendre.jl/issues/130
    ei_broken = (:no_implicit_imports,),
)
