#include<iostream>
#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "RBFMotionSolver.H"
#include "dictionary.H"
#include "fvCFD.H"
#include "IOmanip.H"
#include "SteadyNSSimple.H"
#include <Eigen/Dense>
#define _USE_MATH_DEFINES
#include <cmath>
#include "pointPatchField.H"

class geomPar : public SteadyNSSimple
{
public:
    explicit geomPar(int argc, char* argv[])
        :
        SteadyNSSimple(argc, argv)
    {
        fvMesh& mesh = _mesh();
        IOdictionary* dyndict = new IOdictionary
        (
            IOobject
            (
                "dynamicMeshDictRBF",
                "./constant",
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );
        //ITHACAutilities::getPointsFromPatch(mesh, 6, ahmed, ahmedInd);
        ms = new RBFMotionSolver(mesh, *dyndict);
        vectorField motion(ms->movingPoints().size(), vector::zero);
        movingIDs = ms->movingIDs();
        x0 = ms->movingPoints();
        curX = x0;
        point0 = ms->curPoints();
    }

    List<vector> x0;
    vectorField point0;
    List<vector> curX;
    IOdictionary* dyndict;
    RBFMotionSolver* ms;
    labelList movingIDs;

    labelList getIndicesFromBox(fvMesh& mesh, label ind, Eigen::MatrixXd Box,
                                List<vector>& points2Move)
    {
        points2Move.resize(0);
        pointField meshPoints(mesh.points());
        const polyPatch& patchFound = mesh.boundaryMesh()[ind];
        labelList labelPatchFound(patchFound.meshPoints());
        labelList boxIndices;

        for (int i = 0; i < labelPatchFound.size(); i++)
        {
            auto px = meshPoints[labelPatchFound[i]].component(0);
            auto py = meshPoints[labelPatchFound[i]].component(1);
            auto pz = meshPoints[labelPatchFound[i]].component(2);

            if (px >= min(Box(0, 0), Box(1, 0)) && py >= min(Box(0, 1), Box(1, 1)) &&
                    pz >= min(Box(0, 2), Box(1, 2)) && px <= max(Box(0, 0), Box(1, 0))
                    && py <= max(Box(0, 1), Box(1, 1)) && pz <= max(Box(0, 2), Box(1, 2)) )
            {
                boxIndices.append(labelPatchFound[i]);
                points2Move.append(meshPoints[labelPatchFound[i]]);
            }
        }

        return boxIndices;
    }

    void movePts(double dis, List<vector>& points2Move)
    {
        for (label i = 0; i < points2Move.size(); i++)
        {
            points2Move[i].component(1) = points2Move[i].component(1) + dis;
        }
    }
};

int main(int argc, char* argv[])
{
    geomPar example(argc, argv);
    fvMesh& mesh = example._mesh();
    Eigen::VectorXd beta = Eigen::VectorXd::LinSpaced(50,0,1);
    Eigen::MatrixXd Box(2, 3);

    // Setting the box where you want to serach for boundary
    // Box << x0, y0, z0, x1, y1, z1;
    Box << -7, -7, -1,
        7, 7, 1;
    mkDir("./defGeoms/");
    system("cp -r constant defGeoms/constant");
    system("cp -r system defGeoms/system");
    system("cp -r 0 defGeoms/0");// Box << x0, y0, z0, x1, y1, z1;

    for (int k = 0; k < beta.size(); k++)
    {
        mkDir("./defGeoms/" + name(k + 1));
        mesh.movePoints(example.point0);
        List<vector> points2Move;
        // The cylinder patch is number 5
        labelList boxIndices = example.getIndicesFromBox(mesh, 5, Box, points2Move);
        // Move the points of the cylinder according to a prescribed law (in this case constant vertical displacement)
        example.movePts(beta(k), points2Move);
        ITHACAutilities::setIndices2Value(boxIndices, points2Move, example.movingIDs,
                                          example.curX);
        example.ms->setMotion(example.curX - example.x0);
        // Move any other point (not belonging to the patch) using RBF (see dynamicMeshDictRBF for settings)
        mesh.movePoints(example.ms->curPoints());
        // Write the points in a dedicated folder to visualize geometries
        ITHACAstream::writePoints(mesh.points(), "defGeoms/", name(k + 1) + "/polyMesh/");

    }
}
