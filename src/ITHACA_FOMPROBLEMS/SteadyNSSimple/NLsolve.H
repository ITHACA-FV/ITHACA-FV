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

scalar residual = 1;
scalar uresidual = 1;
Vector<double> uresidual_v(0, 0, 0);

scalar presidual = 1;

scalar csolve = 0;


// Variable that can be changed
turbulence->read();
std::ofstream res_os;
res_os.open("./ITHACAoutput/Offline/residuals", std::ios_base::app);

#if OFVER == 6

while (simple.loop(runTime) && residual > tolerance && csolve < maxIter )
#else
while (simple.loop() && residual > tolerance && csolve < maxIter )
#endif
{
    fvVectorMatrix UEqn(get_Umatrix(U, p));

    if (simple.momentumPredictor())
    {
        Vector<double> uresidual_v = solve(UEqn == -fvc::grad(p)).initialResidual();
        fvOptions.correct(U);
        scalar C = 0;

        for (label i = 0; i < 3; i++)
        {
            if (C < uresidual_v[i])
            {
                C = uresidual_v[i];
            }
        }

        uresidual = C;
    }

    get_Pmatrix(U, p, presidual);
    residual = max(presidual, uresidual);
    Info << "Time = " << runTime.timeName() << nl << endl;
    laminarTransport.correct();
    turbulence->correct();
}

res_os << residual << std::endl;
res_os.close();
runTime.setTime(runTime.startTime(), 0);








