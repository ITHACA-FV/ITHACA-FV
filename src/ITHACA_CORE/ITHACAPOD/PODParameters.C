
#include "ITHACAsystem.H"
#include "PODParameters.H"

template<class Enum>
struct StrLookup: public std::map<std::string,Enum>{
  using Base=std::map<std::string,Enum>;
  using Base::Base;
  Enum lookup(std::string str,Enum defaultValue){
    if (this->count(str)){
      return Base::at(str);
    }else{
      return defaultValue;
    }
  }
};

namespace ITHACAPOD
{

 PODParameters::PODParameters(int argc,char* argv[])
  {

    _args = autoPtr<argList>
        (
          new argList(argc, argv)
          );

    if (!_args->checkRootCase())
    {
      Foam::FatalError.exit();
    }

    argList& args = _args();

    for (int i=0; i<argc; i++)
    {
      Info << "argv[" << i << "] = " << argv[i] << endl;
    }

    runTime0 = autoPtr<Foam::Time>( new Foam::Time( Foam::Time::controlDictName,
                                                    args ) );

    mesh = (
          new fvMesh
          (
            Foam::IOobject
            (
              Foam::fvMesh::defaultRegion,
              runTime0->timeName(),
              *runTime0,
              Foam::IOobject::MUST_READ
              )
            )
          );


    nCells = mesh->cells().size();
    ithacaLibraryParameters = ITHACAparameters::getInstance(*mesh,*runTime0);
    ITHACAdict = ithacaLibraryParameters->ITHACAdict;
    // casenameData = ITHACAdict->lookupOrDefault<fileName>("casename", "./");

    casenameData = ITHACAdict->lookupOrDefault<fileName>("casename", "./");

    fieldlist = static_cast<List<word>>(ITHACAdict->lookup("fields"));
    eigensolver = ithacaLibraryParameters->eigensolver;
    precision = ithacaLibraryParameters->precision;
    exportPython = ithacaLibraryParameters->exportPython;
    exportMatlab = ithacaLibraryParameters->exportMatlab;
    exportTxt = ithacaLibraryParameters->exportTxt;
    outytpe = ithacaLibraryParameters->outytpe;
    pressureResolutionKind = StrLookup<PressureResolutionKind>(
          {
            {"FullOrder",PressureResolutionKind::FullOrder},
            {"ReducedOrder",PressureResolutionKind::ReducedOrder},
            {"Neglected",PressureResolutionKind::Neglected},
          }).lookup(ITHACAdict->lookupOrDefault<word>("pressureResolutionKind", ""),PressureResolutionKind::Undefined);
    nBlocks = ITHACAdict->lookupOrDefault<label>("nBlocks", 1);
    centeredOrNot = ITHACAdict->lookupOrDefault<bool>("centeredOrNot", 1);
    interpFieldCenteredOrNot = ITHACAdict->lookupOrDefault<bool>("interpFieldCenteredOrNot", 0);
    // nMagicPoints = ITHACAdict->lookupOrDefault<label>("nMagicPoints", 1);
    HRMethod = ITHACAdict->lookupOrDefault<word>("HyperReduction", "DEIM");
    HRInterpolatedField = ITHACAdict->lookupOrDefault<word>("DEIMInterpolatedField", "fullStressFunction");
    ECPAlgo = ITHACAdict->lookupOrDefault<word>("ECPAlgo", "Global");
    onLineReconstruct = ITHACAdict->lookupOrDefault<bool>("onLineReconstruct", 0);
    forcingOrNot = ITHACAdict->lookupOrDefault<bool>("forcingOrNot", 0);
    symDiff = ITHACAdict->lookupOrDefault<bool>("symDiff", 0);
    ROMTemporalScheme = ITHACAdict->lookupOrDefault<word>("ROMTemporalScheme", "euler");
           // SOTA can be 0 (no SOTA), D (deterministic version), S (stochastic version)
    useSOTA = ITHACAdict->lookupOrDefault<word>("useSOTA", "None");
    set_useSOTA(useSOTA);

  if (useSOTA != "None")
  {
    Info << "===============================================================" << endl;
    Info << "  RedLUM is launched in " << useSOTA << "-SOTA mode" << endl;
    Info << "  Ithacadict file will be overidden by parameters specified in IthacaFVParameters.C" << endl;
    Info << "===============================================================" << endl;

  }
  
    if ( (!(ROMTemporalScheme == "adams-bashforth"))
         && (!(ROMTemporalScheme == "euler"))
         && (!(ROMTemporalScheme == "euler–maruyama")) )
    {
      Info << "This temporal scheme is not implemented." << endl;
      abort();
    }

    // Object Time to read OpenFOAM data in the correct folder
    runTimeData = new Foam::Time(Foam::Time::controlDictName, ".", casenameData);
    Foam::Time* runTime;
    
    if (Pstream::parRun()) 
    {
      Foam::fileName casenamePar = casenameData + "processor" + name(Pstream::myProcNo());
      Foam::Time* runTimePar = new Foam::Time(Foam::Time::controlDictName, ".", casenamePar);
      runTime = runTimePar;
    }
    else{
      runTime = runTimeData;
    }

    template_field_U = new volVectorField
        (
          IOobject
          (
            "U",
             runTime->times()[1].name(),
            *mesh,
            IOobject::MUST_READ
            ),
          *mesh
          );


    template_field_p = new volScalarField
        (
          IOobject
          (
            "p",
             runTime->times()[1].name(),
            *mesh,
            IOobject::MUST_READ
            ),
          *mesh
          );



    nSimu = ITHACAdict->lookupOrDefault("nSimu", 100);


    // Initialize field_name, field_type and nModes
    field_name.resize(fieldlist.size());
    field_type.resize(fieldlist.size());
    for (label k = 0; k < fieldlist.size(); k++)
    {
      dictionary& subDict = ITHACAdict->subDict(fieldlist[k]);
      field_name[k] = static_cast<word>(subDict.lookup("field_name"));
      field_type[k] = static_cast<word>(subDict.lookup("field_type"));
      nModes.insert(field_name[k],subDict.lookupOrDefault<label>("nmodes",1));
      hilbertSpacePOD.insert(field_name[k],
                             subDict.lookupOrDefault<word>("hilbertSpacePOD","L2"));
      varyingEnergy.insert(field_name[k], 0);
      resolvedVaryingEnergy.insert(field_name[k], 0);
    }

    weightH1 = 1.0;
    weightBC = ITHACAdict->lookupOrDefault<double>("weightBC", 0.);
    patchBC = ITHACAdict->lookupOrDefault<word>("patchBC", "inlet");

    // Initialize startTime, endTime and nSnapshots

    // Get times list from the case folder
    instantList Times = runTime->times();

    // Read Initial and last time from the POD dictionary
    const entry* existnsnap = ITHACAdict->findEntry("Nsnapshots");
    const entry* existLT = ITHACAdict->findEntry("FinalTime");

    // Initiate variable from PODSolverDict
    if ((existnsnap) && (existLT))
    {
      Info << "Error you cannot define LatestTime and NSnapShots together" << endl;
      abort();
    }
    else if (existnsnap)
    {
      InitialTime = ITHACAdict->lookupOrDefault<scalar>("InitialTime", 0);
      FinalTime = ITHACAdict->lookupOrDefault<scalar>("FinalTime", 100000000000000);
      nSnapshots = readScalar(ITHACAdict->lookup("Nsnapshots"));
      startTime = Time::findClosestTimeIndex(runTime->times(), InitialTime);
      nSnapshots = min(nSnapshots , Times.size() - startTime);
      endTime = startTime + nSnapshots - 1;
      FinalTime = std::stof(runTime->times()[endTime].name());
    }
    else
    {
      InitialTime = ITHACAdict->lookupOrDefault<scalar>("InitialTime", 0);
      FinalTime = ITHACAdict->lookupOrDefault<scalar>("FinalTime", 100000000000000);
      endTime = Time::findClosestTimeIndex(runTime->times(), FinalTime);
      startTime = Time::findClosestTimeIndex(runTime->times(), InitialTime);
      nSnapshots = endTime - startTime + 1;
      if (InitialTime > FinalTime)
      {
        Info << "FinalTime cannot be smaller than the InitialTime check your ITHACAdict file\n" << endl;
        abort();
      }
      FinalTime = std::stof(runTime->times()[endTime].name());
    }

    // Read Initial and last time from the POD dictionary
    const entry* existnsnapSimulation = ITHACAdict->findEntry("NsnapshotsSimulation");
    const entry* existLTSimulation = ITHACAdict->findEntry("FinalTimeSimulation");

    scalar InitialTimeSimulation(FinalTime);
    label startTimeSimulation(Time::findClosestTimeIndex(runTime->times(), InitialTimeSimulation));

    if ((existnsnapSimulation) && (existLTSimulation))
    {
      Info << "Error you cannot define LatestTimeSimulation and NSnapShotsSimulation together" << endl;
      abort();
    }
    else if (existnsnapSimulation)
    {
      nSnapshotsSimulation = readScalar(ITHACAdict->lookup("NsnapshotsSimulation"));
      nSnapshotsSimulation = min(nSnapshotsSimulation , Times.size() - startTimeSimulation);
      endTimeSimulation = startTimeSimulation + nSnapshotsSimulation - 1;
      FinalTimeSimulation = std::stof(runTime->times()[endTimeSimulation].name());
    }
    else
    {
      FinalTimeSimulation = ITHACAdict->lookupOrDefault<scalar>("FinalTimeSimulation", 100000000000000);
      endTimeSimulation = Time::findClosestTimeIndex(runTime->times(), FinalTimeSimulation);
      nSnapshotsSimulation = endTimeSimulation - startTimeSimulation + 1;
      if (InitialTimeSimulation > FinalTimeSimulation)
      {
        Info << "FinalTimeSimulation cannot be smaller than the InitialTimeSimulation check your ITHACAdict file\n" << endl;
        abort();
      }
      FinalTimeSimulation = std::stof(runTime->times()[endTimeSimulation].name());
    }

    // Initialize saveTime
    IOdictionary controlDict
        (
          IOobject
          (
            "controlDict",
            runTimeData->system(),
            get_mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
            )
          );
    saveTime = controlDict.lookupOrDefault<double>("writeInterval", 1);

    IOdictionary transportProperties
        (
          IOobject
          (
            "transportProperties",
            runTimeData->constant(),
            *mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
            )
          );

    nu = new dimensionedScalar
        (
          "nu",
          dimViscosity,
          transportProperties
          );

    volume = new volScalarField(
          IOobject(
            "volume",
            runTimeData->timeName(),
            *mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
            ),
          *mesh,
          dimensionedScalar("volume", dimensionSet(0, 3, 0, 0, 0, 0, 0), 1.0)
          );

    volume->primitiveFieldRef() = mesh->V().field();
    Eigen::VectorXd volVect = Foam2Eigen::field2Eigen(mesh->V().field());
    totalVolume = volVect.sum();

    delta = new volScalarField(pow(*volume,1.0/3.0));


    //TO DO : rewrite the following method to search for the object turbulenceProperties and its attibutes (simulationType and LESModel)
    string simulationType;
    string LESModel;

    string DDESModel;

    string filename = "constant/turbulenceProperties";
    std::ifstream strm( filename );

    string line;
    getline( strm, line );

    while (line.find( "simulationType" ) == string::npos || line.find( "//" ) != string::npos)
    {
      getline( strm, line );
    }


    int N = line.size();
    int i =0;
    while ( i<N )
    {
      if (line[i]==' ')
      {
        line.erase(i,1);
        N=N-1;
      }
      else
      {
        i++;
      }
    }

    line.erase(line.size()-1,1);
    line.erase(0,14);
    simulationType=line;



    Info << "--------------------------------------------" << endl;
    Info << "Simulation type:" << simulationType << endl;


    if (simulationType=="laminar")
    {
      Info<< "DNS simulation will be performed" <<endl;
      set_useDNS(true);
      set_useDEIM(false);
    }
    else if(simulationType =="LES")
    {
      template_field_nut = new volScalarField
          (
            IOobject
            (
              "nut",
              runTime->path() + runTime->times()[1].name(),
              *mesh,
              IOobject::MUST_READ
              ),
            *mesh
            );

      getline( strm, line );
      while (line.find( "LESModel" ) == string::npos || line.find( "//" ) != string::npos )
      {
        getline( strm, line );
      }

      N = line.size();
      i =0;
      while ( i<N )
      {
        if (line[i]==' ')
        {
          line.erase(i,1);
          N=N-1;
        }
        else
        {
          i++;
        }
      }
      line.erase(N-1,1);
      line.erase(0,8);
      LESModel=line;


      if (LESModel=="Smagorinsky")
      {
        set_useDEIM(true);
        set_Ck(ITHACAdict->lookupOrDefault<double>("Ck", 0.094));
        set_Ce(ITHACAdict->lookupOrDefault<double>("Ce", 1.048));
        set_useDDES(false);

      }
      else if (LESModel=="kOmegaSSTDDES")
      {
        set_useDDES(true);
        set_useDEIM(false);
        template_field_omega = new volScalarField
            (
              IOobject
              (
                "omega",
                runTime->path() + runTime->times()[1].name(),
                *mesh,
                IOobject::MUST_READ
                ),
              *mesh
              );
        template_field_k = new volScalarField
            (
              IOobject
              (
                "k",
                runTime->path() + runTime->times()[1].name(),
                *mesh,
                IOobject::MUST_READ
                ),
              *mesh
              );
      }
      else
      {
        Info << "Only turbulence model Smagorinsky is supported with simulation type LES" << endl;
        Info << "DDES NOT USED**0**" << endl;


      }


    }
    else{
      Info <<"Simulation type not LES, no DEIM will be performed" << endl;
      set_useDEIM(false);
      Info<<"\n"<<"No DDES was  performed"<<endl;
      set_useDDES(false);
      Info <<"\n"<<"No DNS simulation was perfomed"<<endl;
      set_useDNS(false);
      abort();
    }
    
    if (useDEIM)
    {
      if(get_hilbertSpacePOD()["nut"] == "dL2"){
        set_deltaWeight(ITHACAutilities::getMassMatrixFV(*template_field_nut).array().pow(-2.0/3.0));
      }

      if (HRMethod == "GappyDEIM")
      {
        HRMethod = "DEIM";
      }

      if (interpFieldCenteredOrNot)
      {
        folder_DEIM = "./ITHACAoutput/Hyperreduction/" + HRMethod + "_centered/";
      }
      else
      {
        folder_DEIM = "./ITHACAoutput/Hyperreduction/" + HRMethod + "/";
      }

      if (HRInterpolatedField == "reducedFullStressFunction" || HRInterpolatedField == "reducedNut" )
      {
        word nameToReplace = HRInterpolatedField.substr(7);
        nameToReplace[0] = tolower(nameToReplace[0]);
        HRInterpolatedField += "_" + std::to_string(get_nModes()["U"]) + "modes";
        field_name.append(HRInterpolatedField);
        field_type.append(field_type[find(field_name.begin(), field_name.end(), nameToReplace) - field_name.begin()]);
        nModes.insert(HRInterpolatedField, get_nModes()[nameToReplace]);
        hilbertSpacePOD.insert(HRInterpolatedField,get_hilbertSpacePOD()[nameToReplace]);
        varyingEnergy.insert(HRInterpolatedField, 0);
        resolvedVaryingEnergy.insert(HRInterpolatedField, 0);
      }

      nMagicPoints = get_nModes()[HRInterpolatedField];
      HRSnapshotsField = HRInterpolatedField;
      std::string nbModesUInFolderDEIM = "";
      std::string ECPAlgoInFolderDEIM = "";

      if (HRMethod == "ECP")
      {
        if (!(ECPAlgo == "Global" || ECPAlgo == "EachMode"))
        {
          Info << "Error: ECPAlgo must be Global or EachMode" << endl;
          abort();
        }

        if (HRInterpolatedField == "fullStressFunction")
        {
          HRSnapshotsField = "projFullStressFunction_" + std::to_string(get_nModes()["U"]) + "modes";
        }
        else if (ITHACAutilities::containsSubstring(HRInterpolatedField, "reducedFullStressFunction"))
        {
          HRSnapshotsField = "projReducedFullStressFunction_" + std::to_string(get_nModes()["U"]) + "modes";
        }
        else if (HRInterpolatedField == "nut")
        {
          HRSnapshotsField = "projSmagFromNut_" + std::to_string(get_nModes()["U"]) + "modes";
        }
        else if (ITHACAutilities::containsSubstring(HRInterpolatedField, "reducedNut"))
        {
          HRSnapshotsField = "projSmagFromReducedNut_" + std::to_string(get_nModes()["U"]) + "modes";
        }
        else
        {
          Info << "Error: ECP method is coded only for fullStressFunction and nut" << endl;
          abort();
        }

        nMagicPoints = ITHACAdict->lookupOrDefault("nMagicPoints", nMagicPoints);
        ECPAlgoInFolderDEIM = ECPAlgo + "/";
        nbModesUInFolderDEIM = std::to_string(get_nModes()["U"]) + "modesU/";

        label nECPFields = 1;
        if (ECPAlgo == "EachMode") {nECPFields = nModes["U"];}

        for (label c = 0; c < nECPFields; c++)
        {
          for (label k = 0; k <= (c+1) * (ECPAlgo == "EachMode") * (ITHACAutilities::containsSubstring(HRInterpolatedField,"nut")); k++)
          {
            word fieldNameModec = HRSnapshotsField;
            if (ECPAlgo == "EachMode")
            {
              fieldNameModec += "_" + std::to_string(c+1);
              if (ITHACAutilities::containsSubstring(HRInterpolatedField,"nut")){fieldNameModec += "_" + std::to_string(k);}
            }
            field_name.append(fieldNameModec);
            field_type.append("scalar");
            nModes.insert(fieldNameModec, nModes[HRInterpolatedField]);
            hilbertSpacePOD.insert(field_name.last(),"L2");
            varyingEnergy.insert(field_name.last(), 0);
            resolvedVaryingEnergy.insert(field_name.last(), 0);
          }
        }
      }
      
      folder_DEIM +=  HRInterpolatedField.substr(0,HRInterpolatedField.find("_")) + "/" + nbModesUInFolderDEIM;
      folder_DEIM += std::to_string(nMagicPoints) + "magicPoints/" + ECPAlgoInFolderDEIM;
    }
    else
    {

      folder_DEIM = "./ITHACAoutput/";
    }
  }







Foam::word PODParameters::get_pathHilbertSpace_fromHS(Foam::word hilbertSp)
{
    Foam::word pathHilbertSpace = "";

    if (hilbertSp == "L2" || hilbertSp == "dL2")
    {
        pathHilbertSpace = "";
    }
    else if (hilbertSp == "L2wBC")
    {
        pathHilbertSpace = "_L2wBC";
    }
    else if (hilbertSp == "H1")
    {
        pathHilbertSpace = "_H1";
    }
    else if (hilbertSp == "wH1")
    {
        pathHilbertSpace = "_wH1";
    }
    else
    {
        Foam::Info << "Error: hilbertSpacePOD type " << hilbertSp
                   << " is not valid." << Foam::endl;
        Foam::Info << "dot_product_POD is available for L2, L2wBC, H1 and wH1 only." <<
                   Foam::endl;
        abort();
    }

    return pathHilbertSpace;
}

Foam::word PODParameters::get_pathHilbertSpace(Foam::word fieldName)
{
  return get_pathHilbertSpace_fromHS(hilbertSpacePOD[fieldName]);
}



}
