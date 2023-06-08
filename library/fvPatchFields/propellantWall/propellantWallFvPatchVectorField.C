/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "propellantWallFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "cyclicFvPatch.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::propellantWallFvPatchVectorField::
propellantWallFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    pName_("p"),
    cyclicPatchName_(),
    cyclicPatchLabel_(-1),
    orientation_(1),
    initWallSf_(0),
    initCyclicSf_(0),
    nbrCyclicSf_(0),
    wallFraction_(0),
    curTimeIndex_(-1)
{}


Foam::propellantWallFvPatchVectorField::
propellantWallFvPatchVectorField
(
    const propellantWallFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    pName_(ptf.pName_),
    cyclicPatchName_(ptf.cyclicPatchName_),
    cyclicPatchLabel_(ptf.cyclicPatchLabel_),
    orientation_(ptf.orientation_),
    initWallSf_(ptf.initWallSf_),
    initCyclicSf_(ptf.initCyclicSf_),
    nbrCyclicSf_(ptf.nbrCyclicSf_),
    wallFraction_(ptf.wallFraction_),
    curTimeIndex_(-1)
{}


Foam::propellantWallFvPatchVectorField::
propellantWallFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict, false),
    pName_(dict.getOrDefault<word>("p", "p")),
    cyclicPatchName_(dict.lookup("cyclicPatch")),
    cyclicPatchLabel_(p.patch().boundaryMesh().findPatchID(cyclicPatchName_)),
    orientation_(dict.get<label>("orientation")),
    initWallSf_(p.Sf()),
    initCyclicSf_(p.boundaryMesh()[cyclicPatchLabel_].Sf()),
    nbrCyclicSf_
    (
        refCast<const cyclicFvPatch>
        (
            p.boundaryMesh()[cyclicPatchLabel_],
            dict
        ).neighbFvPatch().Sf()
    ),
    wallFraction_
    (
      initCyclicSf_.size(), 1.0
    ),
    curTimeIndex_(-1)
{
    fvPatchVectorField::operator=(Zero);
}


Foam::propellantWallFvPatchVectorField::
propellantWallFvPatchVectorField
(
    const propellantWallFvPatchVectorField& ptf
)
:
    fixedValueFvPatchVectorField(ptf),
    pName_(ptf.pName_),
    cyclicPatchName_(ptf.cyclicPatchName_),
    cyclicPatchLabel_(ptf.cyclicPatchLabel_),
    orientation_(ptf.orientation_),
    initWallSf_(ptf.initWallSf_),
    initCyclicSf_(ptf.initCyclicSf_),
    nbrCyclicSf_(ptf.nbrCyclicSf_),
    wallFraction_(ptf.wallFraction_),
    curTimeIndex_(-1)
{}


Foam::propellantWallFvPatchVectorField::
propellantWallFvPatchVectorField
(
    const propellantWallFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),
    pName_(ptf.pName_),
    cyclicPatchName_(ptf.cyclicPatchName_),
    cyclicPatchLabel_(ptf.cyclicPatchLabel_),
    orientation_(ptf.orientation_),
    initWallSf_(ptf.initWallSf_),
    initCyclicSf_(ptf.initCyclicSf_),
    nbrCyclicSf_(ptf.nbrCyclicSf_),
    wallFraction_(ptf.wallFraction_),
    curTimeIndex_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::propellantWallFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchVectorField::autoMap(m);

    //- Note: cannot map field from cyclic patch anyway so just recalculate
    //  Areas should be consistent when doing autoMap except in case of
    //  topo changes.
    //- Note: we don't want to use Sf here since triggers rebuilding of
    //  fvMesh::S() which will give problems when mapped (since already
    //  on new mesh)
    const vectorField& areas = patch().boundaryMesh().mesh().faceAreas();
    initWallSf_ = patch().patchSlice(areas);
    initCyclicSf_ = patch().boundaryMesh()
    [
        cyclicPatchLabel_
    ].patchSlice(areas);
    nbrCyclicSf_ = refCast<const cyclicFvPatch>
    (
        patch().boundaryMesh()
        [
            cyclicPatchLabel_
        ]
    ).neighbFvPatch().patch().patchSlice(areas);
}


void Foam::propellantWallFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);

    // See autoMap.
    const vectorField& areas = patch().boundaryMesh().mesh().faceAreas();
    initWallSf_ = patch().patchSlice(areas);
    initCyclicSf_ = patch().boundaryMesh()
    [
        cyclicPatchLabel_
    ].patchSlice(areas);
    nbrCyclicSf_ = refCast<const cyclicFvPatch>
    (
        patch().boundaryMesh()
        [
            cyclicPatchLabel_
        ]
    ).neighbFvPatch().patch().patchSlice(areas);
}


void Foam::propellantWallFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Execute the change to the openFraction only once per time-step
    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        const volScalarField& p = db().lookupObject<volScalarField>
        (
            pName_
        );

        const fvPatch& cyclicPatch = patch().boundaryMesh()[cyclicPatchLabel_];
        const fvPatch& nbrPatch = refCast<const cyclicFvPatch>
        (
            cyclicPatch
        ).neighbFvPatch();

        // Implement OpenFraction algorithm (phi)
        const labelList& fCells = cyclicPatch.faceCells();
        forAll(fCells, cellI)
        {
          if (p[fCells[cellI]] > 0.99) wallFraction_[cellI] = 1.0;
          else wallFraction_[cellI] = p[fCells[cellI]];
        }

        Info<< "wallFraction = " << wallFraction_ << endl;

        vectorField::subField Sfw = this->patch().patch().faceAreas();
        // Info << "wall face cells: " << this->patch().faceCells() << endl;
        // Info << "cyclic face Cells: " << cyclicPatch.faceCells() << endl;
        // Info << "nbr face Cells: " << nbrPatch.faceCells() << endl;
        vectorField newSfw(initWallSf_);
        forAll(newSfw, Sfi)
        {
          if (Sfi <= (fCells.size() - 1))
              newSfw[Sfi] = wallFraction_[Sfi]*newSfw[Sfi];
          else
              newSfw[Sfi] = wallFraction_[Sfi - fCells.size()]*newSfw[Sfi];
        }
        forAll(Sfw, facei)
        {
            Sfw[facei] = newSfw[facei];
        }
        const_cast<scalarField&>(patch().magSf()) = mag(patch().Sf());

        // Update owner side of cyclic
        const_cast<vectorField&>(cyclicPatch.Sf()) =
          max((1.0 - wallFraction_), SMALL)*initCyclicSf_;
        const_cast<scalarField&>(cyclicPatch.magSf()) =
            mag(cyclicPatch.Sf());
        // Update neighbour side of cyclic
        const_cast<vectorField&>(nbrPatch.Sf()) =
            max((1.0 - wallFraction_), SMALL)*nbrCyclicSf_;
        const_cast<scalarField&>(nbrPatch.magSf()) =
            mag(nbrPatch.Sf());

        // Info << "magSf() -> cyclic: " << cyclicPatch.magSf() << endl;
        // Info << "magSf() -> nbr: " << nbrPatch.magSf() << endl;
        // Info << "magSf() -> Wall: " << this->patch().magSf() << endl;

        curTimeIndex_ = this->db().time().timeIndex();
    }

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::propellantWallFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeEntryIfDifferent<word>("p", "p", pName_);
    os.writeEntry("cyclicPatch", cyclicPatchName_);
    os.writeEntry("orientation", orientation_);
    os.writeEntry("openFraction", wallFraction_);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        propellantWallFvPatchVectorField
    );
}


// ************************************************************************* //
