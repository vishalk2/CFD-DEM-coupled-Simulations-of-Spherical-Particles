fEpsilon_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "fEpsilon",
            this->coeffDict_,
            1.0
        )
    ),

    uLim_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "fKupperLimit",
            this->coeffDict_,
            1.0
        )
    ),

    loLim_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "fKlowerLimit",
            this->coeffDict_,
            0.1
        )
    ),

    fK_
    (
        IOobject
        (
            IOobject::groupName("fK", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("zero",loLim_)
    ),

    C2U
    (
        IOobject
        (
            "C2U",
            this->runTime_.timeName(),
            this->mesh_
        ),
        C1_ + (fK_/fEpsilon_)*(C2_ - C1_)
    ),

    delta_
    (
        LESdelta::New
        (
            IOobject::groupName("delta", U.group()),
            *this,
            this->coeffDict_
        )
    ),
    kU_
    (
        IOobject
        (
            IOobject::groupName("kU", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        k_*fK_,
        k_.boundaryField().types()
    ),
    epsilonU_
    (
        IOobject
        (
            IOobject::groupName("epsilonU", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        epsilon_*fEpsilon_,
        epsilon_.boundaryField().types()
    )
{
    bound(k_, this->kMin_);
    bound(epsilon_, this->epsilonMin_);
    bound(kU_, min(fK_)*this->kMin_);
    bound(epsilonU_, fEpsilon_*this->epsilonMin_);


    if (type == typeName)
    {
        this->printCoeffs(type);
        correctNut();
    }
}

    // Update epsilon and G at the wall
    epsilonU_.boundaryFieldRef().updateCoeffs();

    // Unresolved Dissipation equation
    tmp<fvScalarMatrix> epsUEqn
    (
        fvm::ddt(alpha, rho, epsilonU_)
      + fvm::div(alphaRhoPhi, epsilonU_)
      - fvm::laplacian(alpha*rho*DepsilonUEff(), epsilonU_)
     ==
        C1_*alpha()*rho()*G*epsilonU_()/kU_()
      - fvm::SuSp(((2.0/3.0)*C1_ + C3_)*alpha()*rho()*divU, epsilonU_)
      - fvm::Sp(C2U*alpha()*rho()*epsilonU_()/kU_(), epsilonU_)
      + epsilonSource()
      + fvOptions(alpha, rho, epsilonU_)
    );

    epsUEqn.ref().relax();
    fvOptions.constrain(epsUEqn.ref());
    epsUEqn.ref().boundaryManipulate(epsilonU_.boundaryFieldRef());
    solve(epsUEqn);
    fvOptions.correct(epsilonU_);
    bound(epsilonU_, fEpsilon_*this->epsilonMin_);


    // Unresolved Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kUEqn
    (
        fvm::ddt(alpha, rho, kU_)
      + fvm::div(alphaRhoPhi, kU_)
      - fvm::laplacian(alpha*rho*DkUEff(), kU_)
     ==
        alpha()*rho()*G
      - fvm::SuSp((2.0/3.0)*alpha()*rho()*divU, kU_)
      - fvm::Sp(alpha()*rho()*epsilonU_()/kU_(), kU_)
      + kSource()
      + fvOptions(alpha, rho, kU_)
    );

    kUEqn.ref().relax();
    fvOptions.constrain(kUEqn.ref());
    solve(kUEqn);
    fvOptions.correct(kU_);
    bound(kU_, min(fK_)*this->kMin_);

    // Calculation of Turbulent kinetic energy and Dissipation rate
    k_ = kU_/fK_;
    k_.correctBoundaryConditions();
    
    epsilon_ = epsilonU_/fEpsilon_;
    epsilon_.correctBoundaryConditions();
    bound(k_, this->kMin_);
    bound(epsilon_, this->epsilonMin_);

    correctNut();

    // Recalculate fK and C2U with new kU and epsilonU

    // Calculate the turbulence integral length scale
    volScalarField::Internal Lambda
    (
    	pow(k_,1.5)/epsilon_
    );
    
    // update fK
    fK_.primitiveFieldRef() = min(max(
    	sqrt(Cmu_.value())*pow(delta()/Lambda,2.0/3.0), loLim_), uLim_);
    
    // update C2U
    C2U = C1_ + (fK_/fEpsilon_)*(C2_ - C1_);
    
}
