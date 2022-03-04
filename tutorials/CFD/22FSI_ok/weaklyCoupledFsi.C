/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "weaklyCoupledFsi.H"
#include "volFields.H"
#include "Time.H"
#include "IFstream.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
        namespace functionObjects
        {
            defineTypeNameAndDebug(weaklyCoupledFsi, 0);
            
            addToRunTimeSelectionTable(functionObject, weaklyCoupledFsi, dictionary);
        }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::weaklyCoupledFsi::weaklyCoupledFsi
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    forces
    (
        name,
        runTime,
        dict
    ),
    M_(0.0),
    C_(0.0),
    K_(0.0),
    R_(0.0),
    Ymax_(0.0),
    append_(false),
    Y_ (0.0, 0.0),
    Yold_(0.0, 0.0),
    fsiResPtr_(nullptr)
{
    this->read(dict);
    this->createFsiOutFile(dict, append_);
}

Foam::functionObjects::weaklyCoupledFsi::weaklyCoupledFsi
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    forces
    (
        name,
        obr,
        dict
    ),
    M_(0.0),
    C_(0.0),
    K_(0.0),
    R_(0.0),
    Ymax_(0.0),
    append_(false),
    Y_ (0.0, 0.0),
    Yold_(0.0, 0.0),
    fsiResPtr_(nullptr)
{
    this->read(dict);
    createFsiOutFile(dict, append_);
    std::cout << "################ ctor with obr ########################" << std::endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::weaklyCoupledFsi::~weaklyCoupledFsi() {}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::functionObjects::weaklyCoupledFsi::createFsiOutFile(const dictionary& dict, bool appnd)
{
    if (Pstream::master())
    {
        List<word> oldFileLines(0);
        if (appnd)
        {
            IFstream outOld
            (
                dict.get<word>("results")
            );
            
            while(!outOld.eof() && outOld.opened())
            {
                word str(word::null);
                outOld.getLine(str);
                if (!str.empty())
                {
                    oldFileLines.append(str);
                }
            }
            
        }
        
        this->createFsiOutFile(dict);
        
        if (appnd && oldFileLines.size())
        {
            for(label i=1; i<oldFileLines.size(); i++)
            {
                fsiResPtr_() << oldFileLines[i] << endl;
            }
        }
    }
    std::cout << "################ createFsiOutFile method with two variables########################" << std::endl; 
}

void Foam::functionObjects::weaklyCoupledFsi::createFsiOutFile(const dictionary& dict)
{
    if (Pstream::master())
    {
        fsiResPtr_.set
        (
            new OFstream
            (
                dict.get<word>("results")
            )
        );
        fsiResPtr_() << "Time;Y;Vy;Fy" << endl;
    }
    std::cout << "################ createFsiOutFile method with one variable ########################" << std::endl; 

}

bool Foam::functionObjects::weaklyCoupledFsi::read(const dictionary& dict)
{
    if (!forces::read(dict))
    {
        return false;
    }
    
    dict.lookup("M") >> M_;
    
    dict.lookup("C") >> C_;
    
    dict.lookup("K") >> K_;
    
    dict.lookup("R") >> R_;
    
    dict.lookup("Ymax") >> Ymax_;
    
    dict.lookup("append") >> append_;

    Info << "Reading old state" << endl;
    autoPtr<IOdictionary> weaklyCoupledFsiFsiDictPtr;
    //try to read weaklyCoupledFsiFsi object properties
    {
        volVectorField& yDispl = 
            const_cast<volVectorField&>
            (
                obr_.lookupObject<volVectorField>("cellDisplacement") //adapt it with pointDisplacement
            );
            
            //read weaklyCoupledFsiFsiDict header
            IOobject weaklyCoupledFsiFsiHeader
            (
                "weaklyCoupledFsiFsiDict",
                yDispl.mesh().time().timeName(),
                "uniform",
                yDispl.mesh(), // can be seen as objectRegistry
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            );
            
            if (weaklyCoupledFsiFsiHeader.typeHeaderOk<IOdictionary>())
            {
                weaklyCoupledFsiFsiDictPtr.reset
                (
                    new IOdictionary
                    (
                        weaklyCoupledFsiFsiHeader
                    )
                );
                
                Info << "Old state restored" << endl;
                weaklyCoupledFsiFsiDictPtr().lookup("YOld") >> Y_;
                Yold_ = Y_;
            }
            
            setDisplacements(yDispl);
    }
    std::cout << "################ read method call ########################" << std::endl; 
    return true;
}
void Foam::functionObjects::weaklyCoupledFsi::setDisplacements(volVectorField& yDispl)
{
    if (Pstream::parRun())
    {
        Pstream::scatter<scalar>(Y_.first());
    }
    
    vector YPatch (0.0, Y_.first(), 0.0);
    
    forAllConstIter(labelHashSet, patchSet_, iter)
    {
	    label patchId = iter.key();
    	forAll(yDispl.boundaryField()[patchId], faceI)
    	{
    	    yDispl.boundaryFieldRef()[patchId][faceI] = YPatch;
    	}
    }
}

bool Foam::functionObjects::weaklyCoupledFsi::write()
{
    return true;
}

bool Foam::functionObjects::weaklyCoupledFsi::execute()
{
    if (!forces::execute())
    {
        return false;
    }
    
    volVectorField& yDispl = 
        const_cast<volVectorField&>
        (
            obr_.lookupObject<volVectorField>("cellDisplacement")
        );
    
    if (Pstream::master())
    {
        scalar dt = yDispl.mesh().time().deltaT().value();
        scalar ct = yDispl.mesh().time().value();
        
        vector force = forceEff();
        scalar yForce = force.y();
        
        //Runge-Kutta 2-nd order method
        Pair<scalar> Ymid;
        
        Ymid.first() = Yold_.first() + 0.5*dt*Yold_.second();
        Ymid.second()= Yold_.second() + 0.5*dt*(-C_*Yold_.second() - K_*Yold_.first() + R_*yForce) / M_;
        
        Y_.first() = Yold_.first() + dt*Ymid.second();
        Y_.second()= Yold_.second() + dt*(-C_*Ymid.second() - K_*Ymid.first() + R_*yForce) / M_;
        
        if (mag(Y_.first()) >= Ymax_)
        {
            Y_.first() = sign(Y_.first())*Ymax_;
            Y_.second() = (Y_.first() - Yold_.first()) / dt;
        }
        
        Yold_ = Y_;
        
        Log << "yForce = " << yForce << endl;
        Log << "Y= " << Y_.first() << endl;
        Log << "Vy= " << Y_.second() << endl;
        
        fsiResPtr_() << ct << ";" << Y_.first() << ";" << Y_.second() << ";" << yForce << endl;

    }
    
    Pstream::scatter<Pair<scalar> >(Yold_);
    
    //write data to file if time is equal to output time
    if (yDispl.mesh().time().outputTime())
    {
        IOdictionary weaklyCoupledFsiFsiDict
        (
            IOobject
            (
                "weaklyCoupledFsiFsiDict",
                yDispl.mesh().time().timeName(),
                "uniform",
                yDispl.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            )
        );
        
        weaklyCoupledFsiFsiDict.set<Pair<scalar> >
        (
            "YOld",
            Yold_
        );
        
        weaklyCoupledFsiFsiDict.regIOobject::write();
    }
    
    setDisplacements(yDispl);
    
    return true;
}



// ************************************************************************* //
