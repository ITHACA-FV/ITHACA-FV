IOdictionary thermocouplesPosition
(
    IOobject
    (
        "thermocouplesDict",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

List<vector> TCpos(thermocouplesPosition.lookup("positions"));
