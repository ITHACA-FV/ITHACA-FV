import numpy as np
from scipy.linalg import solve
import sys


class ReducedUnsteadyNSExplicit:

    def __init__(self, fluxMethod):
        self.fluxMethod = fluxMethod
        self.N_BC = np.load("./file_python/N_BC.npy")
        self.Nphi_u = np.load("./file_python/Nphi_u.npy")
        self.Nphi_p = np.load("./file_python/Nphi_p.npy")
       
    def solve_online(self, vel):
        if self.fluxMethod == "inconsistent":
            self._solve_inconsistent(vel)
        elif self.fluxMethod == "consistent":
            self._solve_consistent(vel)
        else:
            print("Only the inconsistent and consistent flux methods are implemented.")
            exit(0)

    def _solve_inconsistent(self, vel):

        NUmodes = 10
        NPmodes = 10
        NSUPmodes = 0

        i = 0

        self.tstart = 0
        self.finalTime = 10

        dt = 0.005
        nu = 0.01

        C_tensor = np.load("/home/nrooho/ITHACA-FV/tutorials/CFD/19UnsteadyNSExplicit_turb_caso7_copy/ITHACAoutput/Matrices/C_" + str(0) + "_" + str(NUmodes) + "_" + str(NSUPmodes) + "_t.npy")     

        Cf_tensor = np.load("/home/nrooho/ITHACA-FV/tutorials/CFD/19UnsteadyNSExplicit_turb_caso7_copy/ITHACAoutput/Matrices/Cf_" + str(0) + "_" + str(NUmodes) + "_" + str(NSUPmodes) + "_t.npy")     

        RD_matrix = np.load("/home/nrooho/ITHACA-FV/tutorials/CFD/19UnsteadyNSExplicit_turb_caso7_copy/ITHACAoutput/Matrices/RD/RD" + str(i) + "_" + str(NUmodes) + "_" + str(NSUPmodes) + ".npy")

        RC_matrix = np.load("/home/nrooho/ITHACA-FV/tutorials/CFD/19UnsteadyNSExplicit_turb_caso7_copy/ITHACAoutput/Matrices/RC/RC" + str(i) + "_" + str(NUmodes) + "_" + str(NSUPmodes) + ".npy")

        BP_matrix = np.load("/home/nrooho/ITHACA-FV/tutorials/CFD/19UnsteadyNSExplicit_turb_caso7_copy/ITHACAoutput/Matrices/BP" +  "_" + str(NPmodes) + ".npy")

        P_matrix = np.load("/home/nrooho/ITHACA-FV/tutorials/CFD/19UnsteadyNSExplicit_turb_caso7_copy/ITHACAoutput/Matrices/P" + "_" + str(i) + "_" + str(NUmodes) + "_" + str(NSUPmodes) + "_" + str(NPmodes) + ".npy")

        B_matrix = np.load("/home/nrooho/ITHACA-FV/tutorials/CFD/19UnsteadyNSExplicit_turb_caso7_copy/ITHACAoutput/Matrices/B" + "_" + str(i) + "_" + str(NUmodes) + "_" + str(NSUPmodes) + ".npy")

        K_matrix = np.load("/home/nrooho/ITHACA-FV/tutorials/CFD/19UnsteadyNSExplicit_turb_caso7_copy/ITHACAoutput/Matrices/K" + "_" + str(i) + "_" + str(NUmodes) + "_" + str(NSUPmodes) + "_" + str(NPmodes) + ".npy")


        # LOAD dei file .npy
        a_o = np.load("./file_python/a_o_inc.npy")
        a_n = a_o.copy()   
        b = np.load("./file_python/b_inc.npy")
        x = np.load("./file_python/x_inc.npy")
        presidual = np.load("./file_python/presidual_inc.npy")
        RHS = np.load("./file_python/RHS_inc.npy")

        # Create and resize the solution vectors
        counter = 0
        time = self.tstart
        
        while time < self.finalTime - 0.5 * dt:
            time += dt
            counter += 1
        
        # Set the initial time
        time = self.tstart

        # Set size of online solution
        self.online_solution = [None] * (counter + 1)
        
        # Create vector to store temporal solution and save initial condition as first solution
        a_o = a_o.flatten()
        b = b.flatten()

        tmp_sol = np.zeros((int(self.Nphi_u) + int(self.Nphi_p) + 1))
        tmp_sol[0] = time
        tmp_sol[1: 1 + int(self.Nphi_u)] = a_o          
        tmp_sol[-b.shape[0]:] = b                    
        self.online_solution[0] = tmp_sol

        
        for t in range(1, len(self.online_solution)):
            time += dt
            print(f"################## time = {time} ##################")
            
            M1 = BP_matrix @ a_o * nu
            M2 = P_matrix @ a_o 

            # Pressure Poisson Equation 

            for l in range(int(self.Nphi_p)):
                cf = a_o.T @ Cf_tensor[l, :, :] @ a_o
                RHS[l] = (1 / dt) * M2[l] - cf + M1[l]

            print("########## TERZO SAVE, RIGA 116 ##########")
            print("a_o_in_C", a_o)
            print("b_in_C", b)
            print("x_in_C", x)
            print("presidual_in_C", presidual)
            print("RHS_in_C", RHS)
            print("M1_in_C", M1)
            print("M2_in_C", M2)
            print("cf", cf)

            LinSysDiv = []
            LinSysConv = []
            LinSysDiff = []

            for i in range(int(self.N_BC) + 1):
                
                filename_Div = f"./file_python/LinSysDiv_{i}.npy"
                LinSysDiv.append(np.load(filename_Div))
                
                filename_Conv = f"./file_python/LinSysConv_{i}.npy"
                LinSysConv.append(np.load(filename_Conv))
 
                filename_Diff = f"./file_python/LinSysDiff_{i}.npy"
                LinSysDiff.append(np.load(filename_Diff))

            # Boundary Term (divergence + diffusion + convection)
            RedLinSysP = LinSysDiv.copy()
            RedLinSysP[1] = RHS.copy()

            for i in range(int(self.N_BC)):
                RedLinSysP[1] += vel * ((1 / dt) * LinSysDiv[i + 1] + nu * LinSysDiff[i + 1] +
                                vel * LinSysConv[i + 1])
            
            
            presidual = RedLinSysP[0] @ x - RedLinSysP[1]
            b = np.linalg.solve(RedLinSysP[0], RedLinSysP[1])

            # Momentum Equation

            # Diffusion term 
            M5 = B_matrix @ a_o * nu
            # Pressure gradient term
            M3 = K_matrix @ b

            # boundaryTerm = np.zeros((int(self.Nphi_u), int(self.N_BC)))
            boundaryTerm = np.load("./file_python/boundaryTerm_inc" + ".npy")
            
            for l in range(int(self.N_BC)):
                boundaryTerm[:, l] = vel * (RD_matrix[:, l] * nu + vel * RC_matrix[:, l])

            for l in range(int(self.Nphi_u)):

                cc = a_o.T @ C_tensor[l, :, :] @ a_o
                a_n[l] = a_o[l] + (M5[l] - cc - M3[l]) * dt

                for j in range(int(self.N_BC)):
                    a_n[l] += boundaryTerm[l, j] * dt


            self.online_solution[t] = np.hstack((
                 np.array([time]),       
                 a_n.flatten(),          
                 b.flatten()             
            ))

            print("########## SESTO SAVE, RIGA 207 ##########")
            print("a_o_in_C", a_o)
            print("b_in_C", b)
            print("x_in_C", x)
            print("presidual_in_C", presidual)
            print("RHS_in_C", RHS)
            print("M1_in_C", M1)
            print("M2_in_C", M2)
            print("c_in_C", cc)
            print("a_n_in_C", a_n)

            a_o = a_n.copy()

    def _solve_consistent(self, vel):

        NUmodes = 10
        NPmodes = 10
        NSUPmodes = 0

        i = 0

        self.tstart = 0
        self.finalTime = 10

        dt = 0.005
        nu = 0.01

        C_tensor = np.load("/home/nrooho/ITHACA-FV/tutorials/CFD/19UnsteadyNSExplicit_turb_caso7_copy/ITHACAoutput/Matrices/C_" + str((0)) + "_" + str(NUmodes) + "_" + str(NSUPmodes) + "_t.npy")        

        Cf_tensor = np.load("/home/nrooho/ITHACA-FV/tutorials/CFD/19UnsteadyNSExplicit_turb_caso7_copy/ITHACAoutput/Matrices/Cf_" + str(0) + "_" + str(NUmodes) + "_" + str(NSUPmodes) + "_t.npy")        

        Ci_tensor = np.load("/home/nrooho/ITHACA-FV/tutorials/CFD/19UnsteadyNSExplicit_turb_caso7_copy/ITHACAoutput/Matrices/Ci_" + str(0) + "_" + str(NUmodes) + "_" + str(NSUPmodes) + "_t.npy")        

        RD_matrix = np.load("/home/nrooho/ITHACA-FV/tutorials/CFD/19UnsteadyNSExplicit_turb_caso7_copy/ITHACAoutput/Matrices/RD/RD" + str(i) + "_" + str(NUmodes) + "_" + str(NSUPmodes) + ".npy")

        RC_matrix = np.load("/home/nrooho/ITHACA-FV/tutorials/CFD/19UnsteadyNSExplicit_turb_caso7_copy/ITHACAoutput/Matrices/RC/RC" + str(i) + "_" + str(NUmodes) + "_" + str(NSUPmodes) + ".npy")

        SD_matrix = np.load("/home/nrooho/ITHACA-FV/tutorials/CFD/19UnsteadyNSExplicit_turb_caso7_copy/ITHACAoutput/Matrices/SD/SD" + str(i) + "_" + str(NUmodes) + "_" + str(NSUPmodes) + ".npy")

        SC_matrix = np.load("/home/nrooho/ITHACA-FV/tutorials/CFD/19UnsteadyNSExplicit_turb_caso7_copy/ITHACAoutput/Matrices/SC/SC" + str(i) + "_" + str(NUmodes) + "_" + str(NSUPmodes) + ".npy")

        W_matrix = np.load("/home/nrooho/ITHACA-FV/tutorials/CFD/19UnsteadyNSExplicit_turb_caso7_copy/ITHACAoutput/Matrices/W_" + str(NUmodes) + ".npy")

        BP_matrix = np.load("/home/nrooho/ITHACA-FV/tutorials/CFD/19UnsteadyNSExplicit_turb_caso7_copy/ITHACAoutput/Matrices/BP" +  "_" + str(NPmodes) + ".npy")

        P_matrix = np.load("/home/nrooho/ITHACA-FV/tutorials/CFD/19UnsteadyNSExplicit_turb_caso7_copy/ITHACAoutput/Matrices/P" + "_" + str(i) + "_" + str(NUmodes) + "_" + str(NSUPmodes) + "_" + str(NPmodes) + ".npy")

        B_matrix = np.load("/home/nrooho/ITHACA-FV/tutorials/CFD/19UnsteadyNSExplicit_turb_caso7_copy/ITHACAoutput/Matrices/B" + "_" + str(i) + "_" + str(NUmodes) + "_" + str(NSUPmodes) + ".npy")

        K_matrix = np.load("/home/nrooho/ITHACA-FV/tutorials/CFD/19UnsteadyNSExplicit_turb_caso7_copy/ITHACAoutput/Matrices/K" + "_" + str(i) + "_" + str(NUmodes) + "_" + str(NSUPmodes) + "_" + str(NPmodes) + ".npy")

        I_matrix = np.load("/home/nrooho/ITHACA-FV/tutorials/CFD/19UnsteadyNSExplicit_turb_caso7_copy/ITHACAoutput/Matrices/I" + "_" + str(NUmodes) + ".npy")

        KF_matrix = np.load("/home/nrooho/ITHACA-FV/tutorials/CFD/19UnsteadyNSExplicit_turb_caso7_copy/ITHACAoutput/Matrices/KF" + "_" + str(NUmodes) + "_" + str(NSUPmodes) + ".npy")

        DF_matrix = np.load("/home/nrooho/ITHACA-FV/tutorials/CFD/19UnsteadyNSExplicit_turb_caso7_copy/ITHACAoutput/Matrices/DF" + "_" + str(NUmodes) + "_" + str(NSUPmodes) + ".npy")

        # Create and resize the solution vectors

        a_o = np.load("./file_python/a_o_con" + ".npy")
        a_n = np.zeros(a_o.size)
        b = np.load("./file_python/b_con" + ".npy")
        c_o = np.load("./file_python/c_o_con" + ".npy")
        x = np.load("./file_python/x_con" + ".npy")
        presidual = np.load("./file_python/presidual_con" + ".npy")
        RHS = np.load("./file_python/RHS_con" + ".npy")

        counter = 0
        time = self.tstart
        
        # Number of steps
        while time < self.finalTime - 0.5 * dt:
            time += dt
            counter += 1

        # Set the initial time 
        time = self.tstart

        
        # Set size of online solution 
        self.online_solution = [None] * (counter + 1)
        
        # Create vector to store temporal solution and save initial condition as first solution
        a_o = a_o.flatten()
        b = b.flatten()
        c_o = c_o.flatten()

        tmp_sol = np.zeros((int(self.Nphi_u) + int(self.Nphi_p) + int(self.Nphi_u) + 1))
        tmp_sol[0] = time
        tmp_sol[1:int(self.Nphi_u)+1] = a_o                                                
        tmp_sol[int(self.Nphi_u)+1:int(self.Nphi_u)+1+int(self.Nphi_p)] = b                
        tmp_sol[-int(self.Nphi_u):] = c_o                                               
        self.online_solution[0] = tmp_sol 

        for i in range(1, len(self.online_solution)):
            time += dt
            print(f"################## time = {time} ##################")
            
            # Pressure Poisson Equation

            # Diffusion term 
            M1 = BP_matrix @ a_o * nu
            # Divergence term
            M2 = P_matrix @ a_o 

            for l in range(int(self.Nphi_p)):
                cf = c_o.T @ Cf_tensor[l, :, :] @ a_o
                RHS[l] = (1 / dt) * M2[l] - cf + M1[l]

            print("########## TERZO SAVE, RIGA 116 ##########")
            print("a_o_con_C", a_o)
            print("c_o_con_C", c_o)
            # print("b_con_C", b)
            # print("x_con_C", x)
            # print("presidual_con_C", presidual)
            # print("RHS_con_C", RHS)
            # print("M1_con_C", M1)
            # print("M2_con_C", M2)
            # print("cf_con", cf)

            LinSysDiv = []
            LinSysConv = []
            LinSysDiff = []

            for i in range(int(self.N_BC) + 1):
                
                filename_Div = f"./file_python/LinSysDiv_{i}.npy"
                LinSysDiv.append(np.load(filename_Div))
                
                filename_Conv = f"./file_python/LinSysConv_{i}.npy"
                LinSysConv.append(np.load(filename_Conv))
 
                filename_Diff = f"./file_python/LinSysDiff_{i}.npy"
                LinSysDiff.append(np.load(filename_Diff))

            
            # Boundary Term (divergence + diffusion + convection)
            RedLinSysP = LinSysDiv.copy()
            RedLinSysP[1] = RHS.copy()

            for i in range(int(self.N_BC)):

                RedLinSysP[1] += vel * ((1 / dt) * LinSysDiv[i + 1] + nu * LinSysDiff[i + 1] +
                                   vel * LinSysConv[i + 1])
            

            presidual = RedLinSysP[0] @ x - RedLinSysP[1]
            b = np.linalg.solve(RedLinSysP[0], RedLinSysP[1])

            # Momentum Equation

            # Diffusion term 
            M5 = B_matrix @ a_o * nu
            # Pressure gradient term
            M3 = K_matrix @ b

            print("##########a_o_con_C", a_o)

            boundaryTerm = np.zeros((int(self.Nphi_u), int(self.N_BC)))

            for l in range(int(self.N_BC)):
                boundaryTerm[:, l] = vel * (RD_matrix[:, l] * nu + vel * RC_matrix[:, l])
            
            print("#########a_o_con_C", a_o)
            
            for k in range(int(self.Nphi_u)):
                cc = c_o.T @ C_tensor[k, :, :] @ a_o
                a_n[k] = a_o[k] + (M5[k] - cc - M3[k]) * dt

                for l in range(int(self.N_BC)):
                    a_n[k] += boundaryTerm[k, l] * dt
      
            print("########## SESTO SAVE, RIGA 207 ##########")
            print("a_o_con_C", a_o)
            print("c_o_con_C", c_o)
            # print("a_o_con_C", a_o.shape)
            print("b_con_C", b)
            print("x_con_C", x)
            print("presidual_con_C", presidual)
            print("RHS_con_C", RHS)
            # print("M1_con_C", M1)
            # print("M2_con_C", M2)
            # print("M5_con_C", M5)
            # print("M3_con_C", M3)
            print("cc_con_C", cc)
            print("a_n_con_C", a_n)

            # Flux Equation

            # Mass term
            M6 = I_matrix @ a_o
            # Diffusion term
            M7 = DF_matrix @ a_o * nu
            # Pressure Gradient Term
            M8 = KF_matrix @ b[:, 0]
            # Convective term
            M9 = np.zeros((int(self.Nphi_u)))

            # print("c_o_con_C", c_o)
            # print("c_o_con_C", c_o.shape)

            for k in range(int(self.Nphi_u)):
                M9 += dt * (Ci_tensor[k,:,:] @ a_o).flatten() * c_o[k] 
                # M9 += dt * (Ci_tensor[k,:,:] @ a_o) * c_o[k]

            boundaryTermFlux = np.zeros((int(self.Nphi_u)))

            for l in range(int(self.N_BC)):
                boundaryTermFlux += vel * (SD_matrix[:, l] * nu + vel * SC_matrix[:, l])

            # print("W_matrix", W_matrix.shape)
            # print("boundaryTermFlux", boundaryTermFlux.shape)
            # print("M6", M6.shape)
            # print("M7", M7.shape)
            # print("M8", M8.shape)
            # print("M9", M9.shape)

            # print("W_matrix", W_matrix)
            print("boundaryTermFlux", boundaryTermFlux)
            print("M6", M6)
            print("M7", M7)
            print("M8", M8)
            print("M9", M9)
            
            c_n = np.linalg.solve(W_matrix, M6 - M9 + dt * (-M8 + 
                  M7 + boundaryTermFlux))

            print("c_n", c_n)
            print("c_n.shape", c_n.shape)

            tmp_sol[0] = time 
            tmp_sol[1 : 1 + int(self.Nphi_u)] = a_n.flatten()
            tmp_sol[1 + int(self.Nphi_u):1 + int(self.Nphi_u) + int(self.Nphi_p)] = b.flatten()
            tmp_sol[-int(self.Nphi_u):] = c_n.flatten()
            self.online_solution[i] = tmp_sol

            print("########## SETTIMO SAVE, RIGA 207 ##########")
            print("a_o_con_C", a_o)
            print("b_con_C", b)
            print("x_con_C", x)
            print("presidual_con_C", presidual)
            print("RHS_con_C", RHS)
            # print("M1_con_C", M1)
            # print("M2_con_C", M2)
            # print("M3_con_C", M3)
            # print("M5_con_C", M5)
            # print("M6_con_C", M6)
            # print("M7_con_C", M7)
            # print("M8_con_C", M8)
            # print("M9_con_C", M9)
            print("c_con_C", cc)
            print("c_n", c_n)
            # print("c_n_vec", c_n_vec)
            # print("a_n_con_C", a_n)

            a_o = a_n.copy()
            c_o = c_n

            print("a_o_con_C", a_o)
            print("c_o", c_o)


def reconstruct(self, export_fields, folder):
    if export_fields:
        os.makedirs(folder, exist_ok=True)
        ITHACAutilities.createSymLink(folder)
    print("qui")
    counter = 0
    next_write = 0
    CoeffU = []
    CoeffP = []
    tValues = []
    export_every_index = round(self.exportEvery / self.storeEvery)

    for i in range(len(self.online_solution)):
        if counter == next_write:
            currentUCoeff = self.online_solution[i][1:1 + self.Nphi_u, 0].reshape(-1, 1)
            currentPCoeff = self.online_solution[i][1 + self.Nphi_u:1 + self.Nphi_u + self.Nphi_p, 0].reshape(-1, 1)
            
            CoeffU.append(currentUCoeff)
            CoeffP.append(currentPCoeff)
            
            next_write += export_every_index
            time_now = self.online_solution[i][0, 0]
            tValues.append(time_now)

        counter += 1
    
    uRec = volVectorField("uRec", self.Umodes[0])
    pRec = volScalarField("pRec", self.Pmodes[0])
    
    self.uRecFields = self.Umodes.reconstruct(uRec, CoeffU, "uRec")
    self.pRecFields = self.Pmodes.reconstruct(pRec, CoeffP, "pRec")
    
    if export_fields:
        ITHACAstream.exportFields(self.uRecFields, folder, "uRec")
        ITHACAstream.exportFields(self.pRecFields, folder, "pRec")
