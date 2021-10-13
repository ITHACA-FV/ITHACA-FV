#include "IntegralFilter.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
defineTypeNameAndDebug(IntegralFilter, 0);
addToRunTimeSelectionTable(Filter, IntegralFilter, dictionary);
}

Foam::IntegralFilter::IntegralFilter(const dictionary& dict)
    :
    Filter(),
    dict_(dict)
{
    Info << "Is IntegralFilter" << endl;
}

List<scalarList> Foam::IntegralFilter::apply(const List<labelList>& cells,
        const List<point>& convPoints, const fvMesh& mesh) const
{
    List<scalarList> weights(cells.size());

    for (label i = 0; i < cells.size(); i++)
    {
        scalarList weightsi(cells[i].size());
        scalar factor = 0;

        for (label j = 0; j < cells[i].size(); j++)
        {
            scalar factori = mesh.V()[cells[i][j]] * mag(convPoints[i] - mesh.C()[cells[i][j]]);
            weightsi[j] = factori;
            factor += factori;
        }

        if (factor > VSMALL)
        {
            weights[i] = weightsi / factor;
        }
    }

    return weights;
}
