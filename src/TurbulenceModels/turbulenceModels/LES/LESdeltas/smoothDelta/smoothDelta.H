/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

Class
    Foam::smoothDelta

Description
    Smoothed delta which takes a given simple geometric delta and applies
    smoothing to it such that the ratio of deltas between two cells is no
    larger than a specified amount, typically 1.15.

SourceFiles
    smoothDelta.C

\*---------------------------------------------------------------------------*/

#ifndef smoothDelta_H
#define smoothDelta_H

#include "LESdelta.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                         Class smoothDelta Declaration
\*---------------------------------------------------------------------------*/

class smoothDelta
:
    public LESdelta
{
public:

    //- Public member class used by mesh-wave to propagate the delta-ratio
    class deltaData
    {
        scalar delta_;

        // Private Member Functions

            //- Update. Gets information from neighbouring face/cell and
            //  uses this to update itself (if necessary) and return true.
            template<class TrackingData>
            inline bool update
            (
                const deltaData& w2,
                const scalar scale,
                const scalar tol,
                TrackingData& td
            );


    public:

        // Constructors

            //- Construct null
            inline deltaData();

            //- Construct from delta value
            inline deltaData(const scalar delta);


        // Member Functions

            // Access

            scalar delta() const
            {
                return delta_;
            }


            // Needed by FaceCellWave

                //- Check whether origin has been changed at all or
                //  still contains original (invalid) value.
                template<class TrackingData>
                inline bool valid(TrackingData& td) const;

                //- Check for identical geometrical data.
                //  Used for cyclics checking.
                template<class TrackingData>
                inline bool sameGeometry
                (
                    const polyMesh&,
                    const deltaData&,
                    const scalar,
                    TrackingData& td
                ) const;

                //- Convert any absolute coordinates into relative to
                //  (patch)face centre
                template<class TrackingData>
                inline void leaveDomain
                (
                    const polyMesh&,
                    const polyPatch&,
                    const label patchFacei,
                    const point& faceCentre,
                    TrackingData& td
                );

                //- Reverse of leaveDomain
                template<class TrackingData>
                inline void enterDomain
                (
                    const polyMesh&,
                    const polyPatch&,
                    const label patchFacei,
                    const point& faceCentre,
                    TrackingData& td
                );

                //- Apply rotation matrix to any coordinates
                template<class TrackingData>
                inline void transform
                (
                    const polyMesh&,
                    const tensor&,
                    TrackingData& td
                );

                //- Influence of neighbouring face.
                template<class TrackingData>
                inline bool updateCell
                (
                    const polyMesh&,
                    const label thisCelli,
                    const label neighbourFacei,
                    const deltaData& neighbourInfo,
                    const scalar tol,
                    TrackingData& td
                );

                //- Influence of neighbouring cell.
                template<class TrackingData>
                inline bool updateFace
                (
                    const polyMesh&,
                    const label thisFacei,
                    const label neighbourCelli,
                    const deltaData& neighbourInfo,
                    const scalar tol,
                    TrackingData& td
                );

                //- Influence of different value on same face.
                template<class TrackingData>
                inline bool updateFace
                (
                    const polyMesh&,
                    const label thisFacei,
                    const deltaData& neighbourInfo,
                    const scalar tol,
                    TrackingData& td
                );

                //- Same (like operator==)
                template<class TrackingData>
                inline bool equal(const deltaData&, TrackingData& td) const;

            // Member Operators

                // Needed for List IO
                inline bool operator==(const deltaData&) const;

                inline bool operator!=(const deltaData&) const;

            // IOstream Operators

                friend Ostream& operator<<
                (
                    Ostream& os,
                    const deltaData& wDist
                )
                {
                    return os << wDist.delta_;
                }

                friend Istream& operator>>(Istream& is, deltaData& wDist)
                {
                    return is >> wDist.delta_;
                }
    };


private:

    // Private data

        autoPtr<LESdelta> geometricDelta_;
        scalar maxDeltaRatio_;


    // Private Member Functions

        //- No copy construct and assignment
        smoothDelta(const smoothDelta&) = delete;

        //- No copy assignment
        void operator=(const smoothDelta&) = delete;

        // Calculate the delta values
        void calcDelta();

        //- Fill changedFaces (with face labels) and changedFacesInfo
        //  (with delta).
        //  This is the initial set of faces from which to start the waves.
        //  Since there might be lots of places with delta jumps we can follow
        //  various strategies for this initial 'seed'.
        //  - start from single cell/face and let FaceCellWave pick up all
        //    others from there. might be quite a few waves before everything
        //    settles.
        //  - start from all faces. Lots of initial transfers.
        //  We do something in between:
        //  - start from all faces where there is a jump. Since we cannot easily
        //    determine this across coupled patches (cyclic, processor)
        //    introduce all faces of these and let FaceCellWave sort it out.
        void setChangedFaces
        (
            const polyMesh& mesh,
            const volScalarField& delta,
            DynamicList<label>& changedFaces,
            DynamicList<deltaData>& changedFacesInfo
        );


public:

    //- Runtime type information
    TypeName("smooth");


    // Constructors

        //- Construct from name, turbulenceModel and dictionary
        smoothDelta
        (
            const word& name,
            const turbulenceModel& turbulence,
            const dictionary&
        );


    //- Destructor
    virtual ~smoothDelta()
    {}


    // Member Functions

        //- Read the LESdelta dictionary
        virtual void read(const dictionary&);

        // Correct values
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Data associated with deltaData type are contiguous
template<>
inline bool contiguous<LESModels::smoothDelta::deltaData>()
{
    return true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "smoothDeltaDeltaDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
