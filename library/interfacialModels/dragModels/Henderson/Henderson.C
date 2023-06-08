/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2018 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "Henderson.H"
#include "mathematicalConstants.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace particleDragModels
{
    defineTypeNameAndDebug(Henderson, 0);
    addToRunTimeSelectionTable(particleDragModel, Henderson, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::particleDragModels::Henderson::Henderson
(
    const dictionary& dict,
    const phasePair& pair,
    const bool registerObject
)
:
    particleDragModel(dict, pair, registerObject),
    residualRe_("residualRe", dimless, dict),
    R_
    (
      "R", dimVelocity*dimVelocity/dimTemperature,
      dict.get<scalar>("R")
    ),
    gamma_(dict.get<scalar>("gamma"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::particleDragModels::Henderson::~Henderson()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
Foam::scalar Foam::particleDragModels::Henderson::CdSubsonic
(
  const Foam::scalar& Re,
  const Foam::scalar& M,
  const Foam::scalar& S,
  const Foam::scalar& Tw,
  const Foam::scalar& T
) const
{
  return 24/(Re + S*(4.33 + ((3.65 - 1.53*Tw/T)/(1 + 0.353*Tw/T))*exp(-0.247*Re/S)))
          +  exp(-0.5*M/sqrt(Re))*((4.5 + 0.38*(0.03*Re + 0.48*sqrt(Re)))/(1 + 0.03*Re + 0.48*sqrt(Re)) + 0.1*sqr(M) + 0.2*pow(M, 8))
          + 0.6*S*(1 - exp(-M/Re));
}

Foam::scalar Foam::particleDragModels::Henderson::CdSupersonic
(
  const Foam::scalar& Re,
  const Foam::scalar& M,
  const Foam::scalar& S,
  const Foam::scalar& Tw,
  const Foam::scalar& T
) const
{
  return (0.9 + 0.34/sqr(M) + 1.86*sqrt(M/Re))*(2 + 2/sqr(S) + 1.058*sqrt(Tw/T)/S - 1/pow(S, 4))
          /(1 + 1.86*sqrt(M/Re));
}


Foam::tmp<Foam::volScalarField> Foam::particleDragModels::Henderson::CdRe() const
{
    const tmp<volScalarField> tTw(pair_.dispersed().thermo().T());
    const volScalarField& Tw(tTw());

    const tmp<volScalarField> tT(pair_.continuous().thermo().T());
    const volScalarField& T(tT());

    volScalarField M(max(pair_.magUr()/sqrt(gamma_*R_*T), SMALL));
    volScalarField Re(max(pair_.Re(), SMALL));
    volScalarField S(sqrt(gamma_/2)*M);
    Info << "min -> M: " << min(M).value() << " Re: " << min(Re).value() << " S: " << min(S).value() << endl;
    Info << "max -> M: " << max(M).value() << " Re: " << max(Re).value() << " S: " << max(S).value() << endl;

    volScalarField CdRe
    (
      IOobject
      (
        "CdRe",
        pair_.phase1().mesh()
      ),
      pair_.phase1().mesh(),
      dimensionedScalar("", dimless, 0.0)
    );

    forAll(M, i)
    {
      CdRe[i] = Re[i]*
                (
                  M[i] <= 1.0 ?
                  CdSubsonic(Re[i], M[i], S[i], Tw[i], T[i]) :
                  M[i] >= 1.75 ?
                  CdSupersonic(Re[i], M[i], S[i], Tw[i], T[i]) :
                  (1./3.)*(7. - M[i]*4.)*CdSubsonic(Re[i], 1.0, S[i], Tw[i], T[i])
                  + (4./3.)*(M[i] - 1.0)*(CdSupersonic(Re[i], 1.75, S[i], Tw[i], T[i]))
                );
    }

    CdRe.correctBoundaryConditions();

    // Boundary& bFs(CdRe.boundaryFieldRef());
    // const Boundary& gfbFs = gf.boundaryField();
    // forAll(bFs, i)
    // {
    //   Field<scalar> bFpF = bFs[i];
    //   const Field<scalar> gfbFpF = gfbFs[i];
    //
    //   forAll(bFpF, j)
    //   {
    //     bFpF[j] = std::erf(gfbFpF[j]);
    //   }
    // }

    return Foam::tmp<Foam::volScalarField>(new volScalarField ("tCdRe", CdRe));
}


// ************************************************************************* //
