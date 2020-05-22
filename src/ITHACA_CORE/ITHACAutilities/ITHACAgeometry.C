/*---------------------------------------------------------------------------*\
     ██╗████████╗██╗  ██╗ █████╗  ██████╗ █████╗       ███████╗██╗   ██╗
     ██║╚══██╔══╝██║  ██║██╔══██╗██╔════╝██╔══██╗      ██╔════╝██║   ██║
     ██║   ██║   ███████║███████║██║     ███████║█████╗█████╗  ██║   ██║
     ██║   ██║   ██╔══██║██╔══██║██║     ██╔══██║╚════╝██╔══╝  ╚██╗ ██╔╝
     ██║   ██║   ██║  ██║██║  ██║╚██████╗██║  ██║      ██║      ╚████╔╝
     ╚═╝   ╚═╝   ╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝╚═╝  ╚═╝      ╚═╝       ╚═══╝

 * In real Time Highly Advanced Computational Applications for Finite Volumes
 * Copyright (C) 2017 by the ITHACA-FV authors
-------------------------------------------------------------------------------
License
    This file is part of ITHACA-FV
    ITHACA-FV is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    ITHACA-FV is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU Lesser General Public License for more details.
    You should have received a copy of the GNU Lesser General Public License
    along with ITHACA-FV. If not, see <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/
#include "ITHACAgeometry.H"
using namespace ITHACAutilities;

namespace ITHACAutilities
{

labelList getIndicesFromBox(fvMesh& mesh, List<label> indices,
                            Eigen::MatrixXd Box, List<vector>& points2Move)
{
    points2Move.resize(0);
    labelList boxIndices;
    pointField meshPoints(mesh.points());

    for (label j = 0; j < indices.size(); j++)
    {
        const polyPatch& patchFound = mesh.boundaryMesh()[indices[j]];
        labelList labelPatchFound(patchFound.meshPoints());

        for (int i = 0; i < labelPatchFound.size(); i++)
        {
            auto px = meshPoints[labelPatchFound[i]].component(0);
            auto py = meshPoints[labelPatchFound[i]].component(1);
            auto pz = meshPoints[labelPatchFound[i]].component(2);

            if (px >= min(Box(0, 0), Box(1, 0)) && py >= min(Box(0, 1), Box(1, 1))
                    && pz >= min(Box(0, 2), Box(1, 2)) && px <= max(Box(0, 0), Box(1, 0))
                    && py <= max(Box(0, 1), Box(1, 1)) && pz <= max(Box(0, 2), Box(1, 2)) )
            {
                boxIndices.append(labelPatchFound[i]);
                points2Move.append(meshPoints[labelPatchFound[i]]);
            }
        }
    }

    return boxIndices;
}

List<int> getIndices(fvMesh& mesh, int index, int layers)
{
    List<int> out;
    out.resize(1);
    out[0] = index;

    for (int i = 0; i < layers; i++)
    {
        int size = out.size();

        for (int j = 0; j < size; j++)
        {
            out.append(mesh.cellCells()[out[j]]);
        }
    }

    labelList uniqueIndex;
    uniqueOrder(out, uniqueIndex);
    List<int> out2;
    forAll(uniqueIndex, i)
    {
        out2.append(out[uniqueIndex[i]]);
    }
    return out2;
}

List<int> getIndices(fvMesh& mesh, int index_row,
                     int index_col, int layers)
{
    List<int> out;
    out.resize(2);
    out[0] = index_row;
    out[1] = index_col;

    for (int i = 0; i < layers; i++)
    {
        int size = out.size();

        for (int j = 0; j < size; j++)
        {
            out.append(mesh.cellCells()[out[j]]);
        }
    }

    labelList uniqueIndex;
    uniqueOrder(out, uniqueIndex);
    List<int> out2;
    forAll(uniqueIndex, i)
    {
        out2.append(out[uniqueIndex[i]]);
    }
    return out2;
}

void getPointsFromPatch(fvMesh& mesh, label ind,
                        List<vector>& points, labelList& indices)
{
    pointField meshPoints(mesh.points());
    const polyPatch& patchFound = mesh.boundaryMesh()[ind];
    labelList labelPatchFound(patchFound.meshPoints());
    points.resize(labelPatchFound.size());

    for (int i = 0; i < labelPatchFound.size(); i++)
    {
        points[i] = meshPoints[labelPatchFound[i]];
    }

    indices = labelPatchFound;
}

List<vector> displacedSegment(List<vector> x0, double mux1,
                              double mux2, double muy1, double muy2)
{
    vector minimum = min(x0);
    vector maximum = max(x0);
    double l = Foam::sqrt(pow(minimum[0] - maximum[0],
                              2) + pow(minimum[1] - maximum[1], 2));
    double limin;
    double limax;
    List<vector> xdef(x0);
    List<vector> displ(x0);

    for (int i = 0; i < x0.size(); i++)
    {
        limin = Foam::sqrt(pow(x0[i][0] - minimum[0], 2) + pow(x0[i][1] - minimum[1],
                           2));
        limax = Foam::sqrt(pow(x0[i][0] - maximum[0], 2) + pow(x0[i][1] - maximum[1],
                           2));
        xdef[i][0] = x0[i][0] + mux1 * (1 - limin / l) + mux2 * (1 - limax / l);
        xdef[i][1] = x0[i][1] + muy1 * (1 - limin / l) + muy2 * (1 - limax / l);
        xdef[i][2] = x0[i][2];
    }

    displ = xdef - x0;
    return displ;
}

vector displacePoint(vector x0, vector x_low, vector x_up,
                     double mux_low, double mux_up, double muy_low, double muy_up, double muz_low,
                     double muz_up)
{
    vector direction = x_up - x_low;
    double t0;
    double t1;
    double t2;
    vector x_low_def(x_low[0] + mux_low, x_low[1] + muy_low, x_low[2] + muz_low);
    vector x_up_def(x_up[0] + mux_up, x_up[1] + muy_up, x_up[2] + muz_up);
    vector direction_def = x_up_def - x_low_def;

    if (abs(direction[0]) > 1e-16)
    {
        t0 = (x0[0] - x_low[0]) / direction[0];
    }
    else
    {
        t0 = 0;
    }

    if (abs(direction[1]) > 1e-16)
    {
        t1 = (x0[1] - x_low[1]) / direction[1];

        if (t0 == 0)
        {
            t0 = t1;
        }
    }
    else
    {
        t1 = t0;
    }

    if (abs(direction[2]) > 1e-16)
    {
        t2 = (x0[2] - x_low[2]) / direction[2];

        if (t1 == 0)
        {
            t1 = t2;
        }

        if (t0 == 0)
        {
            t0 = t2;
        }
    }
    else
    {
        t2 = t1;
    }

    Info << abs((x_low + t0 * direction - x0)[0]) << endl;
    Info << abs((x_low + t1 * direction - x0)[1]) << endl;
    Info << abs((x_low + t2 * direction - x0)[2]) << endl;
    M_Assert(abs(abs((x_low + t0 * direction - x0)[0]) + abs((
                     x_low + t0 * direction - x0)[1]) + abs((x_low + t0 * direction - x0)[2])) <
             1e-6, "The givent point is not on the segment");
    vector def_point = x_low_def + t1 * direction_def;
    Info << def_point << endl;
    return def_point;
}

Field<vector> rotateMesh(fvMesh& mesh, double r1, double r2,
                         vector axis,
                         double alpha, labelList movingPointsIDs,
                         List<double> radii, word angleVariationMethod, double v)
{
    M_Assert(angleVariationMethod == "Linear"
             || angleVariationMethod == "Sinusoidal" || angleVariationMethod == "Sigmoid",
             "The variation function of the angle from the inner radius to the outer radius must be either Linear or Sinusoidal or Sigmoid");
    Field<vector> pointRot(mesh.points());

    for (int i = 0; i < movingPointsIDs.size(); i++)
    {
        double l = radii[i];
        vector pointNow = pointRot[movingPointsIDs[i]];

        if (l <= r1)
        {
            double theta = alpha / 180 *  constant::mathematical::pi;
            quaternion q(axis, theta);
            pointRot[movingPointsIDs[i]] = q.transform(pointNow);
        }
        else if (l > r1 && l <= r2)
        {
            if (angleVariationMethod == "Linear")
            {
                double theta = alpha / 180 *  constant::mathematical::pi * (r2 - l) / (r2 - r1);
                quaternion q(axis, theta);
                pointRot[movingPointsIDs[i]] = q.transform(pointNow);
            }
            else if (angleVariationMethod == "Sinusoidal")
            {
                double theta = alpha / 180 *  constant::mathematical::pi * std::sin(
                                   constant::mathematical::pi / 2 * (r2 - l) / (r2 - r1));
                quaternion q(axis, theta);
                pointRot[movingPointsIDs[i]] = q.transform(pointNow);
            }
            else if (angleVariationMethod == "Sigmoid")
            {
                double theta = alpha / 180 *  constant::mathematical::pi * (1 - 1 /
                               (1 + std::exp(
                                    -v * (l - (r1 + r2) / 2))));
                quaternion q(axis, theta);
                pointRot[movingPointsIDs[i]] = q.transform(pointNow);
            }
        }
    }

    return pointRot;
}

Eigen::MatrixXd rotationMatrix(vector AxisOfRotation,
                               double AngleOfRotation)
{
    Eigen::MatrixXd R(3, 3);
    double theta = AngleOfRotation / 180 *  constant::mathematical::pi;
    scalar di = mag(AxisOfRotation);
    scalar ux = AxisOfRotation[0] / di;
    scalar uy = AxisOfRotation[1] / di;
    scalar uz = AxisOfRotation[2] / di;
    R(0, 0) = Foam::cos(theta) + ux * ux * (1 - Foam::cos(theta));
    R(1, 0) = uy * ux * (1 - Foam::cos(theta)) + uz * Foam::sin(theta);
    R(2, 0) = uz * ux * (1 - Foam::cos(theta)) - uy * Foam::sin(theta);
    R(0, 1) = ux * uz * (1 - Foam::cos(theta)) - uz * Foam::sin(theta);
    R(1, 1) = Foam::cos(theta) + uy * uy * (1 - Foam::cos(theta));
    R(2, 1) = ux * uy * (1 - Foam::cos(theta)) + ux * Foam::sin(theta);
    R(0, 2) = ux * uz * (1 - Foam::cos(theta)) + uy * Foam::sin(theta);
    R(1, 2) = uy * uz * (1 - Foam::cos(theta)) - ux * Foam::sin(theta);
    R(2, 2) = Foam::cos(theta) + uz * uz * (1 - Foam::cos(theta));
    return R;
}

Eigen::VectorXd boudaryFaceToCellDistance(
    fvMesh& mesh, label BC_ind)
{
    Eigen::VectorXd cellFaceDistance;
    const polyPatch& cPatch = mesh.boundaryMesh()[BC_ind];
    //Starting index of the face in a patch
    label faceId_start = cPatch.start() ;
    //List of cells close to a boundary
    const labelUList& faceCells = cPatch.faceCells();
    cellFaceDistance.conservativeResize(cPatch.size());
    forAll(cPatch, faceI)
    {
        // index of each face
        label faceID = faceId_start + faceI;
        //id of the owner cell having the face
        label faceOwner = faceCells[faceI] ;
        cellFaceDistance(faceI) = mag(mesh.C()[faceOwner] - mesh.Cf()[faceID]);
    }
    return (cellFaceDistance);
}

List<int> getIndicesFromDisc(fvMesh& mesh, double radius,
                             vector origin, vector axis, List<double>& radii)
{
    pointField meshPoints(mesh.points());
    List<int> indices(meshPoints.size());
    radii.resize(meshPoints.size());
    int k = 0;

    for (int i = 0; i < meshPoints.size(); i++)
    {
        vector pointNow = meshPoints[i];
        vector r = pointNow - origin;
        vector projComponent = (r & axis) / (pow(mag(axis), 2)) * axis;
        vector d = r - projComponent;
        double l = Foam::sqrt(pow(d[0],
                                  2) + pow(d[1], 2) + pow(d[2], 2));

        if (l < radius)
        {
            indices[k] = i;
            radii[k] = l;
            k++;
        }
    }

    indices.resize(k);
    radii.resize(k);
    return indices;
}

template<typename type_f>
List<int> getIndicesFromBox(
    GeometricField<type_f, fvPatchField, volMesh>& field, Eigen::MatrixXd Box)
{
    M_Assert(Box.rows() == 2
             && Box.cols() == 3,
             "The box must be a 2*3 matrix shaped in this way: \nBox = \t|x0, y0, z0|\n\t|x1, yi, z1|\n");
    List<int> indices(field.internalField().size());
    int k = 0;

    for (label i = 0; i < field.internalField().size(); i++)
    {
        auto cx = field.mesh().C()[i][0];
        auto cy = field.mesh().C()[i][1];
        auto cz = field.mesh().C()[i][2];

        if (cx >= Box(0, 0) && cy >= Box(0, 1) && cz >= Box(0, 2) && cx <= Box(1, 0)
                && cy <= Box(1, 1) && cz <= Box(1, 2) )
        {
            indices[k] = i;
            k++;
        }
    }

    indices.resize(k);
    return indices;
}

template List<int> getIndicesFromBox(
    GeometricField<scalar, fvPatchField, volMesh>& field, Eigen::MatrixXd Box);
template List<int> getIndicesFromBox(
    GeometricField<vector, fvPatchField, volMesh>& field, Eigen::MatrixXd Box);

template<typename type_f>
fvMeshSubset* getSubMeshFromBox(
    GeometricField<type_f, fvPatchField, volMesh>& field, Eigen::MatrixXd Box)
{
    List<int> indices = getIndicesFromBox(field, Box);
    fvMeshSubset* sub;
    sub = new fvMeshSubset(field.mesh());
    (field.mesh());
#if OPENFOAM >= 1812
    sub->setCellSubset(indices);
#else
    sub->setLargeCellSubset(indices);
#endif
    return sub;
}

template fvMeshSubset* getSubMeshFromBox(
    GeometricField<scalar, fvPatchField, volMesh>& field, Eigen::MatrixXd Box);
template fvMeshSubset* getSubMeshFromBox(
    GeometricField<vector, fvPatchField, volMesh>& field, Eigen::MatrixXd Box);

volScalarField meshNonOrtho(fvMesh& mesh,
                            volScalarField& NonOrtho)
{
    scalarField sno = (polyMeshTools::faceOrthogonality(mesh, mesh.Sf(), mesh.C()));

    for (int i = 0; i < sno.size(); i++)
    {
        sno[i] = Foam::acos(min(1, sno[i])) * 180 / constant::mathematical::pi;
    }

    surfaceScalarField pippo = mesh.magSf();
    const fvPatchList& patches = mesh.boundary();

    for (int i = 0; i < pippo.internalField().size(); i++)
    {
        pippo.ref()[i] = sno[i];
    }

    for (int i = 0; i < patches.size(); i++)
    {
        if ( patches[i].type() != "empty" )
        {
            label start = patches[i].patch().start();
            label n = patches[i].patch().size();

            for (int k = 0; k < n; k++)
            {
                pippo.boundaryFieldRef()[i][k] = sno[start + k];
            }
        }
    }

    NonOrtho = fvc::average(pippo);
    return NonOrtho;
}

List<vector> rotatePoints(const List<vector>& originalPoints,
                          vector AxisOfRotation, double AngleOfRotation)
{
    double theta = AngleOfRotation / 180 *  constant::mathematical::pi;
    quaternion q(AxisOfRotation, theta);
    List<vector> rotatedPoints(originalPoints);

    for (int i = 0; i < rotatedPoints.size(); i++)
    {
        rotatedPoints[i] = q.transform(rotatedPoints[i]);
    }

    return rotatedPoints;
}
}
