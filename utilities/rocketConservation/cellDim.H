const faceList & ff = mesh.faces();
const pointField & pp = mesh.points();

scalarIOField xDim
(
    IOobject
    (
      "xDim",
      runTime.constant(),
      mesh,
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
    ),
    mesh.C().size()
);

scalarIOField yDim
(
    IOobject
    (
      "yDim",
      runTime.constant(),
      mesh,
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
    ),
    mesh.C().size()
);

scalarIOField zDim
(
    IOobject
    (
      "zDim",
      runTime.constant(),
      mesh,
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
    ),
    mesh.C().size()
);

forAll ( mesh.C(), celli)
{
    const cell & cc = mesh.cells()[celli];
    labelList pLabels(cc.labels(ff));
    pointField pLocal(pLabels.size(), vector::zero);

    forAll (pLabels, pointi)
           pLocal[pointi] = pp[pLabels[pointi]];

    xDim[celli] = Foam::max(pLocal & vector(1,0,0)) - Foam::min(pLocal & vector(1,0,0));
    yDim[celli] = Foam::max(pLocal & vector(0,1,0)) - Foam::min(pLocal & vector(0,1,0));
    zDim[celli] = Foam::max(pLocal & vector(0,0,1)) - Foam::min(pLocal & vector(0,0,1));
}
xDim.write();
yDim.write();
zDim.write();
