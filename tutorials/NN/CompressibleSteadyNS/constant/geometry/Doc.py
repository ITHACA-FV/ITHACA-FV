import numpy as np
import matplotlib.pyplot as plt

# The .txt has to contain a four columns array where the first column is all composed by vs, second column the x values, third the y values and fourth the z values. It gives back a new .txt file with the same structure where new point are there. PAY ATTENTION: the rotation is along the y axe: only x and z are modified. 

f = open('NACA0012_old.obj',"r")
lines=f.readlines()

alpha = -5
alpha = (alpha * np.pi / 180)
#print(np.cos(alpha * np.pi / 180))

x = []
y = []
z = []
y_check = 0

for i in lines:
    if i.split(' ')[0] == "v":
        x.append(i.split(' ')[1])
        y.append(i.split(' ')[2])
        z.append(i.split(' ')[3])
f.close()

for i in range(1,len(y)):
    if abs(float(y[i])) == abs(float(y[i-1])):
        y_check = y_check + abs(float(y[i])) - abs(float(y[i-1]))
        
y_check = int(y_check)

new_x = []
new_y = []
new_z = []
first_vertex = 0

if y_check == 0 :

    print("Rotation around Y axis")

    for i in range(len(x)):
        new_x.append(float(x[i]) * np.cos(alpha) + (-1)*float(z[i]) * np.sin(alpha))
        new_z.append(float(x[i]) * np.sin(alpha) + float(z[i]) * np.cos(alpha))
        
    with open('NACA0012.obj', 'w') as f:
        for i in range(len(lines)):
            if lines[i].split(' ')[0] != "v":
                first_vertex = first_vertex + 1
                f.write(lines[i])
            else :
                f.write('v ' + str(new_x[i-first_vertex]) + ' ' + y[i-first_vertex] + ' ' + str(new_z[i-first_vertex]) + '\n' )
                
else :

    print("Rotation around Z axis")

    for i in range(len(x)):
        new_x.append(float(x[i]) * np.cos(alpha) + (-1)*float(y[i]) * np.sin(alpha))
        new_y.append(float(x[i]) * np.sin(alpha) + float(y[i]) * np.cos(alpha))
        
    with open('NACA0012.obj', 'w') as f:
        for i in range(len(lines)):
            if lines[i].split(' ')[0] != "v":
                first_vertex = first_vertex + 1
                f.write(lines[i])
            else :
                f.write('v ' + str(new_x[i-first_vertex]) + ' ' + new_y[i-first_vertex] + ' ' + str(z[i-first_vertex]) + '\n' ) 

plt.title("New geometry")
plt.axis([0,1,-1,1])
plt.plot(new_x,new_z)
plt.show()

#plt.plot(x,z)
#plt.show()

print ("Done")
