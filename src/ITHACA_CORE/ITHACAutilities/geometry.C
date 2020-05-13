#include "geometry.H"
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

}
