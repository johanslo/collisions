from collisions.dump import *
import matplotlib.pyplot as plt

dump = Dump('./dump.myforce.out')

force, force_coll, positions = dump.find_atom_collsions(425)

dump.list_all_collisions()
# x = []
# y = []
# for i in range(len(dump.timestep)):
#     x.append(dump.data[i][0][0])
#     y.append(dump.data[i][0][1])

# plt.plot(x, y, "k", linestyle = 'dotted')
# markers = "+^s*oD<>v12348pPhHxXd_|"
# for col in dump.collisions:
#     plt.scatter([col.x], [col.y], color = 'k', marker = markers[col.id % len(markers)])

# plt.savefig("3.png")

#print(positions)

c_times = []
c_mag = []
for col in dump.collisions:
    if col.id == 425:
        c_times.append(dump.timestep[col.timestep])
        c_mag.append(force[col.timestep])


plt.xlabel(r'Timestep')
plt.ylabel(r'$f = \sqrt{f_x^2 + f_y^2 + f_z^2}$')
plt.plot(dump.timestep, force, 'k')
plt.savefig("1.png")

plt.close()

plt.xlabel(r'Timestep')
plt.ylabel(r'$f = \sqrt{f_x^2 + f_y^2 + f_z^2}$')
plt.plot(dump.timestep, force_coll, 'k')
plt.scatter(c_times, c_mag, marker='*', color = 'r')
plt.savefig("2.png")

