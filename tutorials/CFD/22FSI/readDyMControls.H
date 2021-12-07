#include "readTimeControls.H"

CorrectPhi_H = pimple.dict().getOrDefault
(
    "correctPhi",
    correctPhi
);

checkMeshCourantNo = pimple.dict().getOrDefault
(
    "checkMeshCourantNo",
    checkMeshCourantNo
);

moveMeshOuterCorrectors = pimple.dict().getOrDefault
(
    "moveMeshOuterCorrectors",
    moveMeshOuterCorrectors
);

//#include "readTimeControls.H"

//correctPhi = pimple.dict().lookupOrDefault("correctPhi", false);

//checkMeshCourantNo = pimple.dict().lookupOrDefault("checkMeshCourantNo", false);
