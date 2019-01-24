#include "ITHACAparallel.H"

ITHACAparallel::ITHACAparallel(fvMesh& mesh)
{
	// Load cell addressing
	indices = new labelIOList(
		IOobject
		(
			"cellProcAddressing",
			mesh.facesInstance(),
			mesh.meshSubDir,
			mesh,
			IOobject::MUST_READ,
			IOobject::NO_WRITE
			)
		); 
	// Load face addressing
	indicesF = new labelIOList(
		IOobject
		(
			"faceProcAddressing",
			mesh.facesInstance(),
			mesh.meshSubDir,
			mesh,
			IOobject::MUST_READ,
			IOobject::NO_WRITE
			)
		); 
	// Calculate total number of cells
	N_IF_glob = (mesh.C().size());
	reduce(N_IF_glob, sumOp<int>()); 

    // BF construction
	Lsize_BF = new labelList (mesh.boundaryMesh().size(),0);
	Gsize_BF = new List<labelList> (mesh.boundaryMesh().size(),labelList(Pstream::nProcs(),0));
	IndFaceLocal = new List<labelList> (mesh.boundaryMesh().size(),labelList(0,0));

	for(int i = 0; i < mesh.boundaryMesh().size(); i++)
	{
		Lsize_BF()[i] = mesh.boundaryMesh()[i].size();
		Gsize_BF()[i][Pstream::myProcNo()] = mesh.boundaryMesh()[i].size();
		reduce(Gsize_BF()[i], sumOp<labelList>());
		IndFaceLocal()[i].resize(Lsize_BF()[i]);
		for (int k = 0; k < mesh.boundaryMesh()[i].size(); k++)
		{
			IndFaceLocal()[i][k] = indicesF()[mesh.boundaryMesh()[i].start()+k];
		}
	}

	Start = new labelList(mesh.boundaryMesh().size(),0);

	for(int i = 0; i < mesh.boundaryMesh().size(); i++)
	{
		if(IndFaceLocal()[i].size() == 0)
		{
			Start()[i] = INT_MAX;	
		}
		else if(IndFaceLocal()[i][0] < 0) 
		{
			Start()[i] = INT_MAX;	
		}
		else
		{
			Start()[i] = IndFaceLocal()[i][0];
		}	
		reduce(Start()[i], minOp<label>()); 
	}
}

void ITHACAparallel::suspendMPI()
{
	Pstream::parRun() = false;
	label comm        = Pstream::worldComm;
	oldProcIDs_       = Pstream::procID(comm);
	newProcIDs_       = List<label> (1, 0);
	Pstream::procID(comm) = newProcIDs_;
}

void ITHACAparallel::resumeMPI()
{
	label comm        = Pstream::worldComm;
	Pstream::procID(comm) = oldProcIDs_;
	Pstream::parRun() = true;
}