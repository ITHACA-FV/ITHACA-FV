Info << "\n Reading nuclear properties \n" << endl;

IOdictionary nuclearProperties
(
    IOobject
    (
        "nuclearProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

_TrefXS = autoPtr<dimensionedScalar>
          (
              new dimensionedScalar
              (
                  nuclearProperties.lookup("TrefXS")
              )
          );
dimensionedScalar& TrefXS = _TrefXS();

_IV1 = autoPtr<dimensionedScalar>
       (
           new dimensionedScalar
           (
               nuclearProperties.lookup("IV1")
           )
       );
dimensionedScalar& IV1 = _IV1();

_D1_0 = autoPtr<dimensionedScalar>
        (
            new dimensionedScalar
            (
                nuclearProperties.lookup("D1_0")
            )
        );
dimensionedScalar& D1_0 = _D1_0();

_alfa_D1 = autoPtr<dimensionedScalar>
           (
               new dimensionedScalar
               (
                   nuclearProperties.lookup("alfa_D1")
               )
           );
dimensionedScalar& alfa_D1 = _alfa_D1();


_A1_0 = autoPtr<dimensionedScalar>
        (
            new dimensionedScalar
            (
                nuclearProperties.lookup("A1_0")
            )
        );
dimensionedScalar& A1_0 = _A1_0();

_alfa_A1 = autoPtr<dimensionedScalar>
           (
               new dimensionedScalar
               (
                   nuclearProperties.lookup("alfa_A1")
               )
           );
dimensionedScalar& alfa_A1 = _alfa_A1();


_NSF1_0 = autoPtr<dimensionedScalar>
          (
              new dimensionedScalar
              (
                  nuclearProperties.lookup("NSF1_0")
              )
          );
dimensionedScalar& NSF1_0 = _NSF1_0();

_alfa_NSF1 = autoPtr<dimensionedScalar>
             (
                 new dimensionedScalar
                 (
                     nuclearProperties.lookup("alfa_NSF1")
                 )
             );
dimensionedScalar& alfa_NSF1 = _alfa_NSF1();

_SP1_0 = autoPtr<dimensionedScalar>
         (
             new dimensionedScalar
             (
                 nuclearProperties.lookup("SP1_0")
             )
         );
dimensionedScalar& SP1_0 = _SP1_0();

_alfa_SP1 = autoPtr<dimensionedScalar>
            (
                new dimensionedScalar
                (
                    nuclearProperties.lookup("alfa_SP1")
                )
            );
dimensionedScalar& alfa_SP1 = _alfa_SP1();

_Keff = autoPtr<dimensionedScalar>
        (
            new dimensionedScalar
            (
                nuclearProperties.lookup("Keff")
            )
        );
dimensionedScalar& Keff = _Keff();

dimensionedScalar K0(Keff);

_K0 = autoPtr<dimensionedScalar>
      (
          new dimensionedScalar(K0)
      );



//read delayed neutron constant

Info << "\n Reading precursors decay constants \n" << endl;

IOdictionary delayedNeutronConstants
(
    IOobject
    (
        "delayedNeutronConstants",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

_lam1 = autoPtr<dimensionedScalar>
        (
            new dimensionedScalar
            (
                delayedNeutronConstants.lookup("lam1")
            )
        );
dimensionedScalar& lam1 = _lam1();

_lam2 = autoPtr<dimensionedScalar>
        (
            new dimensionedScalar
            (
                delayedNeutronConstants.lookup("lam2")
            )
        );
dimensionedScalar& lam2 = _lam2();

_lam3 = autoPtr<dimensionedScalar>
        (
            new dimensionedScalar
            (
                delayedNeutronConstants.lookup("lam3")
            )
        );
dimensionedScalar& lam3 = _lam3();

_lam4 = autoPtr<dimensionedScalar>
        (
            new dimensionedScalar
            (
                delayedNeutronConstants.lookup("lam4")
            )
        );
dimensionedScalar& lam4 = _lam4();

_lam5 = autoPtr<dimensionedScalar>
        (
            new dimensionedScalar
            (
                delayedNeutronConstants.lookup("lam5")
            )
        );
dimensionedScalar& lam5 = _lam5();

_lam6 = autoPtr<dimensionedScalar>
        (
            new dimensionedScalar
            (
                delayedNeutronConstants.lookup("lam6")
            )
        );
dimensionedScalar& lam6 = _lam6();

_lam7 = autoPtr<dimensionedScalar>
        (
            new dimensionedScalar
            (
                delayedNeutronConstants.lookup("lam7")
            )
        );
dimensionedScalar& lam7 = _lam7();

_lam8 = autoPtr<dimensionedScalar>
        (
            new dimensionedScalar
            (
                delayedNeutronConstants.lookup("lam8")
            )
        );
dimensionedScalar& lam8 = _lam8();


_beta1 = autoPtr<dimensionedScalar>
         (
             new dimensionedScalar
             (
                 delayedNeutronConstants.lookup("beta1")
             )
         );
dimensionedScalar& beta1 = _beta1();

_beta2 = autoPtr<dimensionedScalar>
         (
             new dimensionedScalar
             (
                 delayedNeutronConstants.lookup("beta2")
             )
         );
dimensionedScalar& beta2 = _beta2();

_beta3 = autoPtr<dimensionedScalar>
         (
             new dimensionedScalar
             (
                 delayedNeutronConstants.lookup("beta3")
             )
         );
dimensionedScalar& beta3 = _beta3();

_beta4 = autoPtr<dimensionedScalar>
         (
             new dimensionedScalar
             (
                 delayedNeutronConstants.lookup("beta4")
             )
         );
dimensionedScalar& beta4 = _beta4();

_beta5 = autoPtr<dimensionedScalar>
         (
             new dimensionedScalar
             (
                 delayedNeutronConstants.lookup("beta5")
             )
         );
dimensionedScalar& beta5 = _beta5();

_beta6 = autoPtr<dimensionedScalar>
         (
             new dimensionedScalar
             (
                 delayedNeutronConstants.lookup("beta6")
             )
         );
dimensionedScalar& beta6 = _beta6();

_beta7 = autoPtr<dimensionedScalar>
         (
             new dimensionedScalar
             (
                 delayedNeutronConstants.lookup("beta7")
             )
         );
dimensionedScalar& beta7 = _beta7();

_beta8 = autoPtr<dimensionedScalar>
         (
             new dimensionedScalar
             (
                 delayedNeutronConstants.lookup("beta8")
             )
         );
dimensionedScalar& beta8 = _beta8();

