CMAKE_DEPENDENT_OPTION(NEKTAR_SOLVER_HOMOGENEOUS
    "Build the Homogeneous solver." ON
    "NEKTAR_BUILD_SOLVERS" OFF)

IF( NEKTAR_SOLVER_HOMOGENEOUS)
    SET(HomogeneousBloodFlowSolverSource
       ./EquationSystems/CoupledLinearNS.cpp
       ./EquationSystems/CoupledLocalToGlobalC0ContMap.cpp
       ./EquationSystems/IncNavierStokes.cpp
       ./EquationSystems/VelocityCorrectionScheme.cpp
       ./EquationSystems/ViscoElasticConstitutive.cpp
       ./AdvectionTerms/AdjointAdvection.cpp
       ./AdvectionTerms/AdvectionTerm.cpp
       ./AdvectionTerms/LinearisedAdvection.cpp
       ./AdvectionTerms/NavierStokesAdvection.cpp
       ./AdvectionTerms/SkewSymmetricAdvection.cpp
       ./AdvectionTerms/NoAdvection.cpp
       ./ViscoElasticFlowSolver.cpp
       )

    ADD_SOLVER_EXECUTABLE(HomogeneousBloodFlowSolver solvers
        ${HomogeneousBloodFlowSolverSource})


ENDIF( NEKTAR_SOLVER_HOMOGENEOUS )
