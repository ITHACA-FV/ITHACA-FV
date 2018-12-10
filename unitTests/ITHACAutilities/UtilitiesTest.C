#include "ITHACAutilities.H"

bool CreateLink()
{
    ITHACAutilities::createSymLink("./output/pippo");
    return 0;
}

int main(int argc, char **argv)
{
    CreateLink();
    return 0;
}
