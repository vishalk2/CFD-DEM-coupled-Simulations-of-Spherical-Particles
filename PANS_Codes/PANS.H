// Protected data

        // Model coefficients

            dimensionedScalar Cmu_;
            dimensionedScalar C1_;
            dimensionedScalar C2_;
            dimensionedScalar C3_;
            dimensionedScalar sigmak_;
            dimensionedScalar sigmaEps_;
            dimensionedScalar fEpsilon_;
            dimensionedScalar uLim_;
            dimensionedScalar loLim_;

        // Fields

            volScalarField fK_;
            volScalarField C2U;
            //- Run-time selectable delta model
            autoPtr<Foam::LESdelta> delta_;

            volScalarField k_;
            volScalarField kU_;
            volScalarField epsilon_;
            volScalarField epsilonU_;
        
	//- Return the unresolved turbulence kinetic energy
        virtual tmp<volScalarField> kU() const
        {
            return kU_;
        }

        //- Return the unresolved turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilonU() const
        {
            return epsilonU_;
        }
	
	//- Access function to filter width
        inline const volScalarField& delta() const
        {
            return delta_();
        }

        //- Return the effective diffusivity for unresolved k
        tmp<volScalarField> DkUEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DkUEff",
                    (this->nut_/(fK_*fK_*sigmak_/fEpsilon_)
                     + this->nu())
                )
            );
        }

        //- Return the effective diffusivity for unresolved epsilon
        tmp<volScalarField> DepsilonUEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DepsilonUEff",
                    (this->nut_/(fK_*fK_*sigmaEps_/fEpsilon_)
                     + this->nu())
                )
            );
        }
