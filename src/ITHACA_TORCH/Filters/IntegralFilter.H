#ifndef INTEGRALFILTER_H
#define INTEGRALFILTER_H

#include "Filter.H"

namespace Foam
{
class IntegralFilter
    :
    public Filter
{
        //- Disallow copy construct
        IntegralFilter(const IntegralFilter&);

        //- Disallow default bitwise assignment
        void operator=(const IntegralFilter&);

    public:

        //- Runtime type information
        TypeName("IntegralFilter");

        //- Construct from dictionary
        IntegralFilter(const dictionary& dict);

        const dictionary& dict_;


        //- Create and return a clone
        virtual autoPtr<Filter> clone() const
        {
            return autoPtr<Filter>(new IntegralFilter(this->dict_));
        }

        virtual ~IntegralFilter()
        {};

        virtual List<scalarList> apply(const List<labelList>& cells,
                                       const List<point>& convPoints, const fvMesh& mesh) const;


};

}

#endif
