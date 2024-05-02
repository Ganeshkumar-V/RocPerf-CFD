{
  vector Fp = vector(0, 0, 0);
  vector Fs = vector(0, 0, 0);
  scalar Ws = 0;
  scalar Q = 0;

  forAll(mesh.boundary(), bFi)
  {
      // Pressure Force
      const scalarField& pF(p.boundaryField()[bFi]);
      const vectorField& Sf(mesh.boundary()[bFi].Sf());
      // const scalarField& magSf(mesh.boundary()[bFi].magSf());

      Fp += sum(-pF*Sf);

      // Force due to shear stress
      const tensorField& Reffp = tau.boundaryField()[bFi];
      Fs += sum(Sf & Reffp);

      // Workdone due to shear stress
      Ws += sum(Sf & (Reffp & Ugas.boundaryField()[bFi]));

      // Heat Flux
      Q += sum(Sf & qgas.boundaryField()[bFi]);
  }

  // write scalar values to the table
  if (timei%100 == 0)
  {
      file.open(forcefileName, std::ios_base::app);
      file << forceString.str();
      forceString.str(std::string());
      {
        file << runTime.timeName() << ", "
             << Fp.x() << ", "
             << Fp.y() << ", "
             << Fp.z() << ", "
             << mag(Fp) << ", "
             << Fs.x() << ", "
             << Fs.y() << ", "
             << Fs.z() << ", "
             << mag(Fs) << ", "
             << Q << ", "
             << Ws << "\n";
      }
      file.close();
  }
  else
  {
      forceString << runTime.timeName() << ", "
             << Fp.x() << ", "
             << Fp.y() << ", "
             << Fp.z() << ", "
             << mag(Fp) << ", "
             << Fs.x() << ", "
             << Fs.y() << ", "
             << Fs.z() << ", "
             << mag(Fs) << ", "
             << Q << ", "
             << Ws << "\n";
  }
}
