#include "ITHACAsystem.H"
#include "StoredParameters.H"
#include "SnapshotConfiguration.H"

template <class Enum>
struct StrLookup : public std::map<std::string, Enum> {
    using Base = std::map<std::string, Enum>;
    using Base::Base;
    Enum lookup(std::string str, Enum defaultValue)
    {
        if (this->count(str))
        {
            return Base::at(str);
        } else
        {
            return defaultValue;
        }
    }
};


StoredParameters::StoredParameters(int argc, char* argv[]):
    meshConfig_(std::make_unique<MeshConfiguration>())
    , fieldTemplates_(std::make_unique<FieldTemplates>())
    , snapshotConfiguration_(std::make_unique<SnapshotConfiguration>())
    , solverConfiguration_(std::make_unique<SolverConfiguration>())
    , HyperreductionConfiguration_(std::make_unique<HyperreductionConfiguration>())
    , ROMExecutionConfig_(std::make_unique<ROMExecutionConfig>())
    , SimulationFlags_(std::make_unique<SimulationFlags>())
    , PODConfiguration_(std::make_unique<PODConfiguration>())
{

    _args = autoPtr<argList>(
        new argList(argc, argv));

    if (!_args->checkRootCase())
    {
        Foam::FatalError.exit();
    }

    argList& args = _args();

    for (int i = 0; i < argc; i++)
    {
        Info << "argv[" << i << "] = " << argv[i] << endl;
    }

    runTime0 = autoPtr<Foam::Time>(new Foam::Time(Foam::Time::controlDictName,
        args));


    meshConfig_->set_mesh(
        new fvMesh(
            Foam::IOobject(
                Foam::fvMesh::defaultRegion,
                runTime0->timeName(),
                *runTime0,
                Foam::IOobject::MUST_READ)));


    meshConfig_->set_nCells(get_mesh().cells().size());
    ithacaLibraryParameters = ITHACAparameters::getInstance(get_mesh(), *runTime0);
    ITHACAdict = ithacaLibraryParameters->ITHACAdict;

    snapshotConfiguration_->set_casenameData(ITHACAdict->lookupOrDefault<fileName>("casename", "./"));

    PODConfiguration_->set_fieldlist(static_cast<List<word>>(ITHACAdict->lookup("fields")));
    solverConfiguration_->set_eigensolver(ithacaLibraryParameters->eigensolver);
    solverConfiguration_->set_precision(ithacaLibraryParameters->precision);
    SimulationFlags_->set_exportPython(ithacaLibraryParameters->exportPython);
    SimulationFlags_->set_exportMatlab(ithacaLibraryParameters->exportMatlab);
    SimulationFlags_->set_exportTxt(ithacaLibraryParameters->exportTxt);
    solverConfiguration_->set_outytpe(ithacaLibraryParameters->outytpe);
    ROMExecutionConfig_->setPressureResolutionKind(StrLookup<PressureResolutionKind>(
        {
            { "FullOrder", PressureResolutionKind::FullOrder },
            { "ReducedOrder", PressureResolutionKind::ReducedOrder },
            { "Neglected", PressureResolutionKind::Neglected },
        })
            .lookup(ITHACAdict->lookupOrDefault<word>("pressureResolutionKind", ""), PressureResolutionKind::Undefined));
    solverConfiguration_->set_nBlocks(ITHACAdict->lookupOrDefault<label>("nBlocks", 1));
    solverConfiguration_->set_centeredOrNot(ITHACAdict->lookupOrDefault<bool>("centeredOrNot", 1));
    HyperreductionConfiguration_->set_interpFieldCentered(ITHACAdict->lookupOrDefault<bool>("interpFieldCenteredOrNot", 0));
    // nMagicPoints = ITHACAdict->lookupOrDefault<label>("nMagicPoints", 1);
    HyperreductionConfiguration_->setHRMethod(ITHACAdict->lookupOrDefault<word>("HyperReduction", "DEIM"));
    HyperreductionConfiguration_->setHRInterpolatedField(ITHACAdict->lookupOrDefault<word>("DEIMInterpolatedField", "fullStressFunction"));
    HyperreductionConfiguration_->setECPAlgo(ITHACAdict->lookupOrDefault<word>("ECPAlgo", "Global"));
    SimulationFlags_->set_onLineReconstruct(ITHACAdict->lookupOrDefault<bool>("onLineReconstruct", 0));
    SimulationFlags_->set_forcingOrNot(ITHACAdict->lookupOrDefault<bool>("forcingOrNot", 0));
    SimulationFlags_->set_symDiff(ITHACAdict->lookupOrDefault<bool>("symDiff", 0));
    ROMExecutionConfig_->setROMTemporalScheme(ITHACAdict->lookupOrDefault<word>("ROMTemporalScheme", "euler"));
    // SOTA can be 0 (no SOTA), D (deterministic version), S (stochastic version)
    ROMExecutionConfig_->setUseSOTA(ITHACAdict->lookupOrDefault<word>("useSOTA", "None"));

    if (ROMExecutionConfig_->useSOTA() != "None")
    {
        Info << "===============================================================" << endl;
        Info << "  RedLUM is launched in " << ROMExecutionConfig_->useSOTA() << "-SOTA mode" << endl;
        Info << "  Ithacadict file will be overidden by parameters specified in m_parameters.C" << endl;
        Info << "===============================================================" << endl;
    }

    if ((!(ROMExecutionConfig_->ROMTemporalScheme() == "adams-bashforth")) && (!(ROMExecutionConfig_->ROMTemporalScheme() == "euler")) && (!(ROMExecutionConfig_->ROMTemporalScheme() == "euler–maruyama")))
    {
        Info << "This temporal scheme is not implemented." << endl;
        abort();
    }

    // Object Time to read OpenFOAM data in the correct folder
    runTimeData = new Foam::Time(Foam::Time::controlDictName, ".", get_casenameData());
    Foam::Time* runTime;

    if (Pstream::parRun())
    {
        Foam::fileName casenamePar = get_casenameData() + "processor" + name(Pstream::myProcNo());
        Foam::Time* runTimePar = new Foam::Time(Foam::Time::controlDictName, ".", casenamePar);
        runTime = runTimePar;
    } else
    {
        runTime = runTimeData;
    }


    fieldTemplates_->set_U(new volVectorField(
        IOobject(
            "U",
            runTime->times()[1].name(),
            get_mesh(),
            IOobject::MUST_READ),
        get_mesh()));


    fieldTemplates_->set_p(new volScalarField(
        IOobject(
            "p",
            runTime->times()[1].name(),
            get_mesh(),
            IOobject::MUST_READ),
        get_mesh()));


    snapshotConfiguration_->set_nSimu(ITHACAdict->lookupOrDefault("nSimu", 100));


    initializeFieldConfiguration();

    PODConfiguration_->set_weightH1(1.0);
    PODConfiguration_->set_weightBC(ITHACAdict->lookupOrDefault<double>("weightBC", 0.));
    PODConfiguration_->set_patchBC(ITHACAdict->lookupOrDefault<word>("patchBC", "inlet"));

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
    } else if (existnsnap)
    {
        snapshotConfiguration_->set_initialTime(ITHACAdict->lookupOrDefault<scalar>("InitialTime", 0));
        snapshotConfiguration_->set_finalTime(ITHACAdict->lookupOrDefault<scalar>("FinalTime", 100000000000000));
        snapshotConfiguration_->set_nSnapshots(readScalar(ITHACAdict->lookup("Nsnapshots")));
        snapshotConfiguration_->set_startTime(Time::findClosestTimeIndex(runTime->times(), get_initialTime()));
        snapshotConfiguration_->set_nSnapshots(min(get_nSnapshots(), Times.size() - get_startTime()));
        snapshotConfiguration_->set_endTime(get_startTime() + get_nSnapshots() - 1);
        snapshotConfiguration_->set_finalTime(std::stof(runTime->times()[get_endTime()].name()));
    } else
    {
        snapshotConfiguration_->set_initialTime(ITHACAdict->lookupOrDefault<scalar>("InitialTime", 0));
        snapshotConfiguration_->set_finalTime(ITHACAdict->lookupOrDefault<scalar>("FinalTime", 100000000000000));
        snapshotConfiguration_->set_endTime(Time::findClosestTimeIndex(runTime->times(), get_finalTime()));
        snapshotConfiguration_->set_startTime(Time::findClosestTimeIndex(runTime->times(), get_initialTime()));
        snapshotConfiguration_->set_nSnapshots(get_endTime() - get_startTime() + 1);
        if (get_initialTime() > get_finalTime())
        {
            Info << "FinalTime cannot be smaller than the InitialTime check your ITHACAdict file\n"
                 << endl;
            abort();
        }
        snapshotConfiguration_->set_finalTime(std::stof(runTime->times()[get_endTime()].name()));
    }

    // Read Initial and last time from the POD dictionary
    const entry* existnsnapSimulation = ITHACAdict->findEntry("NsnapshotsSimulation");
    const entry* existLTSimulation = ITHACAdict->findEntry("FinalTimeSimulation");

    scalar InitialTimeSimulation(get_finalTime());
    label startTimeSimulation(Time::findClosestTimeIndex(runTime->times(), InitialTimeSimulation));

    if ((existnsnapSimulation) && (existLTSimulation))
    {
        Info << "Error you cannot define LatestTimeSimulation and NSnapShotsSimulation together" << endl;
        abort();
    } else if (existnsnapSimulation)
    {
        snapshotConfiguration_->set_nSnapshotsSimulation(readScalar(ITHACAdict->lookup("NsnapshotsSimulation")));
        snapshotConfiguration_->set_nSnapshotsSimulation(min(get_nSnapshotsSimulation(), Times.size() - startTimeSimulation));
        snapshotConfiguration_->set_endTimeSimulation(startTimeSimulation + get_nSnapshotsSimulation() - 1);
        snapshotConfiguration_->set_finalTimeSimulation(std::stof(runTime->times()[get_endTimeSimulation()].name()));
    } else
    {
        snapshotConfiguration_->set_finalTimeSimulation(ITHACAdict->lookupOrDefault<scalar>("FinalTimeSimulation", 100000000000000));
        snapshotConfiguration_->set_endTimeSimulation(Time::findClosestTimeIndex(runTime->times(), get_finalTimeSimulation()));
        snapshotConfiguration_->set_nSnapshotsSimulation(get_endTimeSimulation() - startTimeSimulation + 1);
        if (InitialTimeSimulation > get_finalTimeSimulation())
        {
            Info << "FinalTimeSimulation cannot be smaller than the InitialTimeSimulation check your ITHACAdict file\n"
                 << endl;
            abort();
        }
        snapshotConfiguration_->set_finalTimeSimulation(std::stof(runTime->times()[get_endTimeSimulation()].name()));
    }

    // Initialize writeInterval
    IOdictionary controlDict(
        IOobject(
            "controlDict",
            runTimeData->system(),
            get_mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE));
    snapshotConfiguration_->set_writeInterval(controlDict.lookupOrDefault<double>("writeInterval", 1));

    IOdictionary transportProperties(
        IOobject(
            "transportProperties",
            runTimeData->constant(),
            get_mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE));

    ROMExecutionConfig_->setNu(new dimensionedScalar(
        "nu",
        dimViscosity,
        transportProperties));

    meshConfig_->set_volume(new volScalarField(
        IOobject(
            "volume",
            runTimeData->timeName(),
            get_mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE),
        get_mesh(),
        dimensionedScalar("volume", dimensionSet(0, 3, 0, 0, 0, 0, 0), 1.0)));

    meshConfig_->get_volume().primitiveFieldRef() = get_mesh().V().field();
    Eigen::VectorXd volVect = Foam2Eigen::field2Eigen(get_mesh().V().field());
    meshConfig_->set_totalVolume(volVect.sum());

    meshConfig_->set_delta(new volScalarField(pow(get_volume(), 1.0 / 3.0)));


    // TO DO : rewrite the following method to search for the object turbulenceProperties and its attibutes (simulationType and LESModel)
    string simulationType;
    string LESModel;

    string DDESModel;

    string filename = "constant/turbulenceProperties";
    std::ifstream strm(filename);

    string line;
    getline(strm, line);

    while (line.find("simulationType") == string::npos || line.find("//") != string::npos)
    {
        getline(strm, line);
    }


    int N = line.size();
    int i = 0;
    while (i < N)
    {
        if (line[i] == ' ')
        {
            line.erase(i, 1);
            N = N - 1;
        } else
        {
            i++;
        }
    }

    line.erase(line.size() - 1, 1);
    line.erase(0, 14);
    simulationType = line;


    Info << "--------------------------------------------" << endl;
    Info << "Simulation type:" << simulationType << endl;


    if (simulationType == "laminar")
    {
        Info << "DNS simulation will be performed" << endl;
        set_useDNS(true);
        set_useDEIM(false);
    } else if (simulationType == "LES")
    {
        fieldTemplates_->set_nut(new volScalarField(
            IOobject(
                "nut",
                runTime->path() + runTime->times()[1].name(),
                get_mesh(),
                IOobject::MUST_READ),
            get_mesh()));

        getline(strm, line);
        while (line.find("LESModel") == string::npos || line.find("//") != string::npos)
        {
            getline(strm, line);
        }

        N = line.size();
        i = 0;
        while (i < N)
        {
            if (line[i] == ' ')
            {
                line.erase(i, 1);
                N = N - 1;
            } else
            {
                i++;
            }
        }
        line.erase(N - 1, 1);
        line.erase(0, 8);
        LESModel = line;


        if (LESModel == "Smagorinsky")
        {
            set_useDEIM(true);
            set_Ck(ITHACAdict->lookupOrDefault<double>("Ck", 0.094));
            set_Ce(ITHACAdict->lookupOrDefault<double>("Ce", 1.048));
            set_useDDES(false);

        } else if (LESModel == "kOmegaSSTDDES")
        {
            set_useDDES(true);
            set_useDEIM(false);
            fieldTemplates_->set_omega(new volScalarField(
                IOobject(
                    "omega",
                    runTime->path() + runTime->times()[1].name(),
                    get_mesh(),
                    IOobject::MUST_READ),
                get_mesh()));
            fieldTemplates_->set_k(new volScalarField(
                IOobject(
                    "k",
                    runTime->path() + runTime->times()[1].name(),
                    get_mesh(),
                    IOobject::MUST_READ),
                get_mesh()));
        } else
        {
            Info << "Only turbulence model Smagorinsky is supported with simulation type LES" << endl;
            Info << "DDES NOT USED**0**" << endl;
        }


    } else
    {
        Info << "Simulation type not LES, no DEIM will be performed" << endl;
        set_useDEIM(false);
        Info << "\n"
             << "No DDES was  performed" << endl;
        set_useDDES(false);
        Info << "\n"
             << "No DNS simulation was perfomed" << endl;
        set_useDNS(false);
        abort();
    }

    HyperreductionConfiguration_->setFolderDEIM("./ITHACAoutput/");
    if (SimulationFlags_->useDEIM())
    {
        if (get_hilbertSpacePOD()["nut"] == "dL2")
        {
            set_deltaWeight(ITHACAutilities::getMassMatrixFV(fieldTemplates_->get_nut()).array().pow(-2.0 / 3.0));
        }
        HyperreductionConfiguration_->initializeHyperreduction(PODConfiguration_);
    }
}

void StoredParameters::initializeFieldConfiguration()
{
    const auto& fieldlist = PODConfiguration_->fieldlist();
    // Resize field name/type arrays
    Foam::List<Foam::word> field_names(fieldlist.size());
    Foam::List<Foam::word> field_types(fieldlist.size());
    // Read from dictionary and populate PODModeConfig
    for (label k = 0; k < fieldlist.size(); k++)
    {
        dictionary& subDict = ITHACAdict->subDict(fieldlist[k]);

        field_names[k] = static_cast<word>(subDict.lookup("field_name"));
        field_types[k] = static_cast<word>(subDict.lookup("field_type"));

        label nMode = subDict.lookupOrDefault<label>("nmodes", 1);
        word hilbertSpace = subDict.lookupOrDefault<word>("hilbertSpacePOD", "L2");

        // Populate the PODModeConfig module
        PODConfiguration_->insert_nModes(field_names[k], nMode);
        PODConfiguration_->insert_hilbertSpacePOD(field_names[k], hilbertSpace);
        PODConfiguration_->set_varyingEnergy(field_names[k], 0);
        PODConfiguration_->set_resolvedVaryingEnergy(field_names[k], 0);
    }

    // Store field names and types
    PODConfiguration_->set_field_name(field_names);
    PODConfiguration_->set_field_type(field_types);
}


Foam::word StoredParameters::get_pathHilbertSpace_fromHS(Foam::word hilbertSp)
{
    Foam::word pathHilbertSpace = "";

    if (hilbertSp == "L2" || hilbertSp == "dL2")
    {
        pathHilbertSpace = "";
    } else if (hilbertSp == "L2wBC")
    {
        pathHilbertSpace = "_L2wBC";
    } else if (hilbertSp == "H1")
    {
        pathHilbertSpace = "_H1";
    } else if (hilbertSp == "wH1")
    {
        pathHilbertSpace = "_wH1";
    } else
    {
        Foam::Info << "Error: hilbertSpacePOD type " << hilbertSp
                   << " is not valid." << Foam::endl;
        Foam::Info << "dot_product_POD is available for L2, L2wBC, H1 and wH1 only." << Foam::endl;
        abort();
    }

    return pathHilbertSpace;
}

Foam::word StoredParameters::get_pathHilbertSpace(Foam::word fieldName)
{
    return get_pathHilbertSpace_fromHS(PODConfiguration_->hilbertSpacePOD()[fieldName]);
}

