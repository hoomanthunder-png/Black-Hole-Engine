sun = 1989000000000000000000000000000
c = 299792458
c2 = 89875517873681764
c3 = 26944002417373989539335912
gravitational_constant = 0.000000000066743
pi = 3.1415926
time_type_list = {"s" : "Second", "m" : "Minute", "h" : "Hour", "d" : "Day", "mn" : "Month", "y" : "Year"}
distance_type_list = {"m" : 1, "km" : 1000, "au" : 149597870700, "ly" : 9460730000000000}
velocity_type_list = {"cm/s" : 0.01, "m/s" : 1, "sos/s" : 343, "km/s" : 1000, "c/s" : 299792458}
length_type_list = [
    (1, lambda l: l * 1000, "MilliMeter"),
    (1000, lambda l: l, "Meter"),
    (149597870700, lambda l: l / 1000, "KiloMeter"),
    (9460730000000000, lambda l: l / 149597870700, "Astronomical Unit"),
    (float("inf"), lambda l: l / 9460730000000000, "Light Year"),
]

speed_type_list = [
    (1, lambda s: s * 100, "cm/s"),
    (100, lambda s: s, "m/s"),
    (1000, lambda s: s / 343, "sos/s"),
    (30000, lambda s: s / 1000, "km/s"),
    (float("inf"), lambda s: s / 299792458, "c")
]

color_hex_list = [
    (0, "1"),
    (1, "1"),
    (2, "2"),
    (3, "3"),
    (4, "4"),
    (5, "5"),
    (6, "6"),
    (7, "7"),
    (8, "8"),
    (9, "9"),
    (10, "A"),
    (11, "B"),
    (12, "C"),
    (13, "D"),
    (14, "E"),
    (15, "F"),
]

from math import sqrt, cos, acos, radians, sin
import NumberLex
import warnings
import numpy as np
from colorama import Fore
import scipy.optimize as opt

#matplotlib liberary
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.patheffects as pe
fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")
fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
plt.show(block=False)

warnings.filterwarnings("ignore", category=RuntimeWarning)

print(Fore.MAGENTA + "Black Hole Engine. Version = 14.719" + Fore.RESET)
print()

def main():
    exit_first_loop = 0
    time_type = ""
    distance_type = ""
    #getting black hole type
    while True:
        bh_type = input(Fore.LIGHTYELLOW_EX + "Black Hole Type (schw or kerr):\n"
            + Fore.LIGHTBLACK_EX + "schw or 1 = non-Spinning Black Hole\n"
            "kerr or 2 = Spinning Black Hole\n"
            + Fore.RESET + "==============================\n"
            ">>>").strip()
        if type_check(bh_type) == False:
            continue
        elif type_check(bh_type) == "schw":
            a_star = 0
            print(Fore.LIGHTGREEN_EX + "Type: " + Fore.RESET + "Non-Spinning Black Hole", end="\n")
            print()
            break
        elif type_check(bh_type) == "kerr":
            print(Fore.LIGHTGREEN_EX + "Type: " + Fore.RESET + "Spinning Black Hole", end="\n")
            print()
            break

    #getting black hole's spin parameter
    if bh_type == "kerr" or bh_type == "2":
        while True:
            a_star = input(Fore.LIGHTYELLOW_EX + "What is the Dimensionless Spin Parameter (must be between 0 and 1):\n"
                + Fore.LIGHTBLACK_EX + "0 = The Black Hole is Actually a Schwarzschild Black Hole.\n"
                "1 = The Black Hole has the Maximum Theoretical Spin Speed!\n"
                + Fore.RESET + "==============================\n"
                ">>>").strip()

            result = spin_parser(a_star)

            if result[0] == False:
                continue
            elif result[0] == None:
                bh_type = "schw"
                a_star = 0
                print(Fore.LIGHTGREEN_EX + "Type: " + Fore.RESET + "Non-Spinning Black Hole", end="\n")
                print()
                break
            else:
                a_star = result[0]
                velocity = result[1]

            print(Fore.LIGHTGREEN_EX + f"Dimensionless Spin Parameter: {Fore.RESET} {a_star}", end="\n")
            print(Fore.LIGHTGREEN_EX + f"Black Hole Velocity at the Equator: {Fore.RESET} {speed_unit_manager(velocity)}", end="\n")
            print()
            break

    #getting black hole mass
    while True:
        bh_mass = input(Fore.LIGHTYELLOW_EX + "Mass Of The Black Hole:\n"
            + Fore.LIGHTBLACK_EX + "Enter a Number in order to Sun Masses\n"
            "or Add 'kg' at the end\n"
            + Fore.RESET + "==============================\n"
            ">>>").strip().lower()

        mass_result = mass_parser(bh_mass)
        if mass_result[0] == False:
            continue
        else:
            solar_mass = mass_result[0]
            bh_mass = mass_result[1]
            black_hole_info_print(bh_type, bh_mass, solar_mass, a_star, False, False, False, False, False)
            break

    #select, print, and construct custome informations
    bh_info_list = select_bh_info(bh_type)
    black_hole_info_print(bh_type, bh_mass, solar_mass, a_star, bh_info_list[0], bh_info_list[1], bh_info_list[2], bh_info_list[3], bh_info_list[4])
    black_hole_construct(bh_type, bh_mass, a_star, bh_info_list[0], bh_info_list[1], bh_info_list[2], bh_info_list[3], bh_info_list[4])

    #select, and construct space-time
    st_type_list = select_space_time_type(bh_type)
    space_time_construct(bh_type, bh_mass, a_star, st_type_list[0], st_type_list[1], st_type_list[2], st_type_list[3], outer_event_horizon_radius(bh_type, bh_mass, a_star, np.radians(90), False))

    #select, print, and construct photon orbitals
    orbital_list = select_photon_orbitals(bh_type)
    photon_orbitals_info_print(bh_type, bh_mass, a_star, orbital_list[0], orbital_list[1], orbital_list[2], orbital_list[3], orbital_list[4], orbital_list[5], orbital_list[6], orbital_list[7], orbital_list[8], orbital_list[9])
    photon_orbitals_construct(bh_type, bh_mass, a_star, orbital_list[0], orbital_list[1], orbital_list[2], orbital_list[3], orbital_list[4], orbital_list[5], orbital_list[6], orbital_list[7], orbital_list[8], orbital_list[9])

    #select, print, and construct ISCO orbitals
    orbital_list = select_isco_orbitals(bh_type)
    isco_orbitals_info_print(bh_type, bh_mass, a_star, orbital_list[0], orbital_list[1], orbital_list[2], orbital_list[3], orbital_list[4], orbital_list[5], orbital_list[6], orbital_list[7], orbital_list[8], orbital_list[9])
    isco_orbitals_construct(bh_type, bh_mass, a_star, orbital_list[0], orbital_list[1], orbital_list[2], orbital_list[3], orbital_list[4], orbital_list[5], orbital_list[6], orbital_list[7], orbital_list[8], orbital_list[9])

    #select, print, and construct metter/light experience
    xp_list = select_object_experience_type(bh_type, bh_mass, a_star)
    object_experience(bh_type, bh_mass, a_star, xp_list[0], xp_list[1], xp_list[2], xp_list[3], xp_list[4], xp_list[5], xp_list[6], xp_list[7])

    input("Press Enter to exit...")

    # #getting type of the experience
    # while True:
    #     near_type = input(Fore.LIGHTYELLOW_EX + "Alright! you want to be Stationary or Orbiting?\n"
    #     + Fore.RESET + "==============================\n"
    #     ">>>").strip().lower()
    #     match near_type:
    #         case "stationary":
    #             while True:
    #                 print()
    #                 time = input(Fore.LIGHTYELLOW_EX + "how long you want to be near it?\n"
    #                 + Fore.LIGHTBLACK_EX + "use 's' as Seconds\n"
    #                 "use 'm' as Minutes\n"
    #                 "use 'h' as Hours\n"
    #                 "use 'd' as Days\n"
    #                 "use 'mn' as Months\n"
    #                 "use 'y' as Years\n"
    #                 "or enter without any to use 'h' as default\n"
    #                 + Fore.RESET + "==============================\n"
    #                 ">>>")
    #                 for i in time:
    #                     if i.isalpha():
    #                         time_type = time_type + str(i) 
    #                         time = time.replace(i, "")

    #                 if time_type == "": time_type = "h"
    #                 time = int(time)
    #                 dilated_time = time_dilation(a_star, "stationary", "", bh_mass, distance, time, 90)
    #                 if dilated_time == False:
    #                     print(Fore.RED + "You're In the Black Hole's Ergosphere. You Can't be Stationary in here!!!" + Fore.RESET)
    #                 else:
    #                     print(Fore.LIGHTGREEN_EX + f"If you spend {Fore.LIGHTBLUE_EX} {time:,.3f} {time_type_list[time_type]} {Fore.LIGHTGREEN_EX} near this Black Hole while you're Stationary, then about {Fore.LIGHTBLUE_EX} {dilated_time:,.3f} {time_type_list[time_type]} {Fore.LIGHTGREEN_EX} will pass in the universe!", end="\n")
    #                     print()
    #                     print(Fore.LIGHTCYAN_EX + f"so you'll Time Travel about {Fore.LIGHTBLUE_EX} {dilated_time - time:,.3f} {time_type_list[time_type]}" + Fore.RESET)             
    #                 exit_first_loop = 1
    #                 break
    #         case "orbiting":
    #             if bh_type == "schw":
    #                 #calculating orbital parameters
    #                 print()
    #                 print(Fore.RED + "as a Rotating Object around a Non-Spinning Black Hole:" + Fore.RESET)
    #                 print("========================================")
    #                 if r_isco(a_star, "", gravitational_radius(bh_mass), "", False) <= distance:
    #                     print(Fore.LIGHTBLUE_EX + "Safe Distance and Stable Orbit!" + Fore.RESET)
    #                     print()
    #                 elif photon_sphere_radius(a_star, gravitational_radius(bh_mass), False) <= distance < r_isco(a_star, "prograde", gravitational_radius(bh_mass), "", False):
    #                     print(Fore.RED + "your'e between Photon Sphere and ISCO so: " + Fore.RESET)
    #                     print(Fore.LIGHTGREEN_EX + "you must have: ")
    #                     print()
    #                 else:
    #                     print(Fore.RED + "your'e Inside the Photon Sphere so: " + Fore.RESET)
    #                     print(Fore.LIGHTGREEN_EX + "you must have: ")
    #                     print()
                    
    #                 orbital_angular_velocity(bh_mass, 0, gravitational_radius(bh_mass), distance, "", True)
    #                 orbital_period(orbital_angular_velocity(bh_mass, 0, gravitational_radius(bh_mass), distance, "", False), True)
    #                 orbital_frequensy(orbital_angular_velocity(bh_mass, 0, gravitational_radius(bh_mass), distance, "", False), True)
    #                 coordinate_orbital_velocity(orbital_angular_velocity(bh_mass, 0, gravitational_radius(bh_mass), distance, "", False), distance, True)
    #                 local_orbital_velocity(a_star, event_horizon_radius(bh_type, bh_mass, 0, False)[1], orbital_angular_velocity(bh_mass, 0, gravitational_radius(bh_mass), distance, "",False), frame_dragging_angular_velocity(bh_mass, 0, distance, gravitational_radius(bh_mass), False), distance, gravitational_radius(bh_mass), True)

    #                 if photon_sphere_radius(a_star, gravitational_radius(bh_mass), False) <= distance < r_isco(a_star, "prograde", gravitational_radius(bh_mass), "", False):
    #                     print(Fore.LIGHTRED_EX + "which is Nearly Impossible!" + Fore.RESET)
    #                 elif distance < photon_sphere_radius(a_star, gravitational_radius(bh_mass), False):
    #                     print(Fore.LIGHTRED_EX + "then it's Absolutely Impossible!" + Fore.RESET)
    #                 print()

    #                 while True:                    
    #                     time = input(Fore.LIGHTYELLOW_EX + "how long you want to orbit near it?\n"
    #                     + Fore.LIGHTBLACK_EX + "use 's' as Seconds\n"
    #                     "use 'm' as Minutes\n"
    #                     "use 'h' as Hours\n"
    #                     "use 'd' as Days\n"
    #                     "use 'mn' as Months\n"
    #                     "use 'y' as Years\n"
    #                     "or enter without any to use 'h' as default\n"
    #                     + Fore.RESET + "==============================\n"
    #                     ">>>")
    #                     for i in time:
    #                         if i.isalpha():
    #                             time_type = time_type + str(i) 
    #                             time = time.replace(i, "")

    #                     if time_type == "": time_type = "h"
    #                     time = int(time)
    #                     dilated_time = time_dilation(a_star, "orbiting", "", bh_mass, distance, time, "")
    #                     if dilated_time == False:
    #                         print(Fore.RED + "There is no stable orbit at this distance!" + Fore.RESET)
    #                     else:
    #                         print(Fore.LIGHTGREEN_EX + f"If you spend {Fore.LIGHTBLUE_EX} {time:,.3f} {time_type_list[time_type]} {Fore.LIGHTGREEN_EX} near this Black Hole while you're Orbiting, then about {Fore.LIGHTBLUE_EX} {dilated_time:,.3f} {time_type_list[time_type]} {Fore.LIGHTGREEN_EX} will pass in the universe!", end="\n")
    #                         print()
    #                         print(Fore.LIGHTCYAN_EX + f"so you'll Time Travel about {Fore.LIGHTBLUE_EX} {dilated_time - time:,.3f} {time_type_list[time_type]}"+ Fore.RESET)
    #                     exit_first_loop = 1
    #                     break
    #             else:
    #                 while True:
    #                     print()
    #                     direction = input(Fore.LIGHTYELLOW_EX + "Prograde or Retrograde? " + Fore.RESET).strip().lower()
    #                     match direction:
    #                         case "prograde":
    #                             #calculating prograde orbital parameters
    #                             print()
    #                             print(Fore.RED + "as a Co_Rotating Object around a Spinning Black Hole:" + Fore.RESET)
    #                             print("========================================")
    #                             if r_isco(a_star, "prograde", gravitational_radius(bh_mass), "", False) <= distance and a_star != 1:
    #                                 print(Fore.LIGHTBLUE_EX + "Safe Distance and Stable Orbit!" + Fore.RESET)
    #                                 print()
    #                             elif photon_sphere_radius(a_star, gravitational_radius(bh_mass), False)[0] < distance < r_isco(a_star, "prograde", gravitational_radius(bh_mass), "", False):
    #                                 print(Fore.RED + "your'e between Co-Rotating Photon Sphere and ISCO so: " + Fore.RESET)
    #                                 print(Fore.LIGHTGREEN_EX + "you must have: ")
    #                                 print()
    #                             else:
    #                                 print(Fore.RED + "your'e Inside the Co-Rotating Photon Sphere so: " + Fore.RESET)
    #                                 print(Fore.LIGHTGREEN_EX + "you must have: ")
    #                                 print()
                                    
    #                             orbital_angular_velocity(bh_mass, a_star, gravitational_radius(bh_mass), distance, "prograde", True)
    #                             orbital_period(orbital_angular_velocity(bh_mass, a_star, gravitational_radius(bh_mass), distance, "prograde", False), True)
    #                             orbital_frequensy(orbital_angular_velocity(bh_mass, a_star, gravitational_radius(bh_mass), distance, "prograde", False), True)
    #                             coordinate_orbital_velocity(orbital_angular_velocity(bh_mass, a_star, gravitational_radius(bh_mass), distance, "prograde", False), distance, True)
    #                             local_orbital_velocity(a_star, event_horizon_radius(bh_type, bh_mass, a_star, False)[1], orbital_angular_velocity(bh_mass, a_star, gravitational_radius(bh_mass), distance, "prograde",False), frame_dragging_angular_velocity(bh_mass, a_star, distance, gravitational_radius(bh_mass), False), distance, gravitational_radius(bh_mass), True)

    #                             if photon_sphere_radius(a_star, gravitational_radius(bh_mass), False)[0] <= distance < r_isco(a_star, "prograde", gravitational_radius(bh_mass), "", False):
    #                                 print(Fore.LIGHTRED_EX + "which is Nearly Impossible!" + Fore.RESET)
    #                             elif distance <= photon_sphere_radius(a_star, gravitational_radius(bh_mass), False)[0]:
    #                                 print(Fore.LIGHTRED_EX + "then it's Absolutely Impossible!" + Fore.RESET)
    #                             print()

    #                             while True:                    
    #                                 time = input(Fore.LIGHTYELLOW_EX + "how long you want to orbit prograde near it?\n"
    #                                 + Fore.LIGHTBLACK_EX + "use 's' as Seconds\n"
    #                                 "use 'm' as Minutes\n"
    #                                 "use 'h' as Hours\n"
    #                                 "use 'd' as Days\n"
    #                                 "use 'mn' as Months\n"
    #                                 "use 'y' as Years\n"
    #                                 "or enter without any to use 'h' as default\n"
    #                                 + Fore.RESET + "==============================\n"
    #                                 ">>>")
    #                                 #calculating time dilation for prograde orbit
    #                                 for i in time:
    #                                     if i.isalpha():
    #                                         time_type = time_type + str(i)
    #                                         time = time.replace(i, "")

    #                                 if time_type == "": time_type = "h"
    #                                 time = int(time)
    #                                 dilated_time = time_dilation(a_star, "orbiting", "prograde", bh_mass, distance, time, 90)
    #                                 if dilated_time == False:
    #                                     print(Fore.RED + "There is no stable orbit at this distance!" + Fore.RESET)
    #                                 else:
    #                                     print(Fore.LIGHTGREEN_EX + f"If you spend {Fore.LIGHTBLUE_EX} {time:,.3f} {time_type_list[time_type]} {Fore.LIGHTGREEN_EX} near this Black Hole while you're Orbiting Prograde, then about {Fore.LIGHTBLUE_EX} {dilated_time:,.3f} {time_type_list[time_type]} {Fore.LIGHTGREEN_EX} will pass in the universe!", end="\n")
    #                                     print(Fore.LIGHTCYAN_EX + f"so you'll Time Travel about {Fore.LIGHTBLUE_EX} {dilated_time - time:,.3f} {time_type_list[time_type]}"+ Fore.RESET)
    #                                 exit_first_loop = 1
    #                                 break
    #                         case "retrograde":
    #                             #calculating retrograde orbital parameters
    #                             print()
    #                             print(Fore.RED + "as a Counter_Rotating Object around a Spinning Black Hole:" + Fore.RESET)
    #                             print("========================================")
    #                             if r_isco(a_star, "retrograde", gravitational_radius(bh_mass), "", False) <= distance:
    #                                 print(Fore.LIGHTBLUE_EX + "Safe Distance and Stable Orbit!" + Fore.RESET)
    #                                 print()
    #                             elif photon_sphere_radius(a_star, gravitational_radius(bh_mass), False)[1] < distance < r_isco(a_star, "retrograde", gravitational_radius(bh_mass), "", False):
    #                                 print(Fore.RED + "your'e between Counter-Rotating Photon Sphere and ISCO so: " + Fore.RESET)
    #                                 print(Fore.LIGHTGREEN_EX + "you must have: ")
    #                                 print()
    #                             else:
    #                                 print(Fore.RED + "your'e Inside the Counter-Rotating Photon Sphere so: " + Fore.RESET)
    #                                 print(Fore.LIGHTGREEN_EX + "you must have: ")
    #                                 print()

    #                             orbital_angular_velocity(bh_mass, a_star, gravitational_radius(bh_mass), distance, "retrograde", True)
    #                             orbital_period(orbital_angular_velocity(bh_mass, a_star, gravitational_radius(bh_mass), distance, "retrograde", False), True)
    #                             orbital_frequensy(orbital_angular_velocity(bh_mass, a_star, gravitational_radius(bh_mass), distance, "retrograde", False), True)
    #                             coordinate_orbital_velocity(orbital_angular_velocity(bh_mass, a_star, gravitational_radius(bh_mass), distance, "retrograde", False), distance, True)
    #                             local_orbital_velocity(a_star, event_horizon_radius(bh_type, bh_mass, a_star, False)[1], orbital_angular_velocity(bh_mass, a_star, gravitational_radius(bh_mass), distance, "retrograde",False), frame_dragging_angular_velocity(bh_mass, a_star, distance, gravitational_radius(bh_mass), False), distance, gravitational_radius(bh_mass), True)

                                
    #                             if photon_sphere_radius(a_star, gravitational_radius(bh_mass), False)[1] <= distance < r_isco(a_star, "retrograde", gravitational_radius(bh_mass), "", False):
    #                                 print(Fore.LIGHTRED_EX + "which is Nearly Impossible!" + Fore.RESET)
    #                             elif distance < photon_sphere_radius(a_star, gravitational_radius(bh_mass), False)[1]:
    #                                 print(Fore.LIGHTRED_EX + "then it's Absolutely Impossible!" + Fore.RESET)
    #                             print()

    #                             while True:                    
    #                                 time = input(Fore.LIGHTYELLOW_EX + "how long you want to orbit retrograde near it?\n"
    #                                 + Fore.LIGHTBLACK_EX + "use 's' as Seconds\n"
    #                                 "use 'm' as Minutes\n"
    #                                 "use 'h' as Hours\n"
    #                                 "use 'd' as Days\n"
    #                                 "use 'mn' as Months\n"
    #                                 "use 'y' as Years\n"
    #                                 "or enter without any to use 'h' as default\n"
    #                                 + Fore.RESET + "==============================\n"
    #                                 ">>>")
    #                                 #calculating time dilation for retrograde orbit
    #                                 for i in time:
    #                                     if i.isalpha():
    #                                         time_type = time_type + str(i)
    #                                         time = time.replace(i, "")

    #                                 if time_type == "": time_type = "h"
    #                                 time = int(time)
    #                                 dilated_time = time_dilation(a_star, "orbiting", "retrograde", bh_mass, distance, time, 90)
    #                                 if dilated_time == False:
    #                                     print(Fore.RED + "There is no stable orbit at this distance!"+ Fore.RESET)
    #                                 else:
    #                                     print(Fore.LIGHTGREEN_EX + f"If you spend {Fore.LIGHTBLUE_EX} {time:,.3f} {time_type_list[time_type]} {Fore.LIGHTGREEN_EX} near this Black Hole while you're Orbiting Retrograde, then about {Fore.LIGHTBLUE_EX} {dilated_time:,.3f} {time_type_list[time_type]} {Fore.LIGHTGREEN_EX} will pass in the universe!", end="\n")
    #                                     print()
    #                                     print(Fore.LIGHTCYAN_EX + f"so you'll Time Travel about {Fore.LIGHTBLUE_EX} {dilated_time - time:,.3f} {time_type_list[time_type]}"+ Fore.RESET)
    #                                 exit_first_loop = 1
    #                                 break
    #                         case _:
    #                             print(Fore.RED + "must be prograde or retrograde!" + Fore.RESET)
    #                             continue
    #                     if exit_first_loop == 1:
    #                         break
    #         case _:
    #             print(Fore.RED + "must be stationary or orbiting" + Fore.RESET)
    #             continue
    #     if exit_first_loop == 1:
    #         break

#getting and checking the black hole's 3 first important parameters
def type_check(bh_type):
    bh_type.lower()
    bh_type.strip()

    if not bh_type in ["schw", "1", "kerr", "2"]:
        print(f"{Fore.LIGHTCYAN_EX} type_check says:{Fore.LIGHTRED_EX} Incorrect bh_type! must be 'schw | 1' or 'kerr | 2' !")
        print()
        return False
    elif bh_type in ["schw", "1"]:
        return "schw"
    elif bh_type in ["kerr", "2"]:
        return "kerr"

def spin_parser(a_star):
    a_star.lower()
    unit = ""
    velocity = 0

    if a_star == "":
        a_star = "0"

    if a_star == "0":
        print(Fore.LIGHTGREEN_EX + "by setting the Spin Parameter as Zero, you actually chose a Non_Spinning Black Hole!" + Fore.RESET)
        print()
        return None, None

    if a_star.startswith("."):
        a_star = "0" + a_star

    first_dot = a_star.find(".")

    if first_dot != -1 and "c" not in a_star:
        a_star = (a_star[:first_dot + 1] + a_star[first_dot + 1:].replace(".", ""))

        part1, part2 = a_star.split(".")

        for i in part1:
            if not i.isdigit():
                part1 = part1.replace(str(i), "")

        for i in part2:
            if not i.isdigit():
                part2 = part2.replace(str(i), "")

        a_star = part1 + "." + part2
        a_star = float(a_star)

        if 0 <= a_star <= 1:
            velocity = c * (1 - np.sqrt(1 - a_star**2)) / a_star
            return a_star, velocity
        else:
            print(
                f"{Fore.LIGHTCYAN_EX} spin parser says:"
                f"{Fore.LIGHTRED_EX} Incorrect spin! dimensionless spin must be between 0 and 1!"
            )
            print()
            return False, None

    elif first_dot != -1 and "c" in a_star:
        a_star = (a_star[:first_dot + 1] + a_star[first_dot + 1:].replace(".", ""))

        first_c = a_star.find("c")

        if first_c != -1:
            a_star = a_star[:first_c]

        for i in a_star:
            if not i.isdigit() and i != ".":
                a_star = a_star.replace(str(i), "")

        a_star = float(a_star)

        if 0 <= a_star <= 1:
            velocity = a_star * c
            a_star = (2 * a_star) / (1 + a_star**2)
            return a_star, velocity
        else:
            print(
                f"{Fore.LIGHTCYAN_EX} spin parser says:"
                f"{Fore.LIGHTRED_EX} Incorrect spin! speed must be less or equal to c!"
            )
            print()
            return False, None
        
    elif first_dot == -1 and "c" in a_star:
        for i in a_star:
            if i.isalpha():
                a_star = a_star.replace(str(i), "")

        if a_star != "":
            a_star = float(a_star)
        else:
            a_star = 1
        
        if 0 <= a_star <= 1:
            velocity = a_star * c
            a_star = (2 * a_star) / (1 + a_star**2)
            return a_star, velocity
        else:
            print(
                f"{Fore.LIGHTCYAN_EX} spin parser says:"
                f"{Fore.LIGHTRED_EX} Incorrect spin! speed must be less or equal to c!"
            )
            print()
            return False, None
        
    else:
        for i in a_star:
            if not i.isdigit():
                unit += str(i)

        if unit == "" and 0 <= float(a_star) <= 1:
            a_star = float(a_star)
            velocity = c * (1 - np.sqrt(1 - a_star**2)) / a_star
            return a_star, velocity

        elif unit == "" and not 0 <= float(a_star) <= 1:
            print(
                f"{Fore.LIGHTCYAN_EX} spin parser says:"
                f"{Fore.LIGHTRED_EX} Incorrect spin! dimensionless spin must be between 0 and 1!"
            )
            print()
            return False, None

        elif unit != "" and not unit in ["km/s", "kms", "m/s", "ms", "c"]:
            print(
                f"{Fore.LIGHTCYAN_EX} spin parser says:"
                f"{Fore.LIGHTRED_EX} Incorrect spin! unit must be either 'km/s' | 'kms' | 'm/s' | 'ms' | 'c'!"
            )
            print()
            return False, None

        else:
            a_star = float(a_star)

            if unit in ["km/s", "kms"]:
                velocity = a_star * 1000
            elif unit in ["m/s", "ms"]:
                velocity = a_star
            elif unit == "c":
                velocity = a_star * c

            if 0 <= velocity <= c:
                a_star = (2 * (velocity / c)) / (1 + (velocity / c)**2)
                return a_star, velocity
            else:
                print(
                    f"{Fore.LIGHTCYAN_EX} spin parser says:"
                    f"{Fore.LIGHTRED_EX} Incorrect spin! speed must be less or equal to "
                    f"299,792,458 m/s | 299,792.458 km/s | c!"
                )
                print()
                return False, None

def mass_parser(bh_mass):
    if bh_mass != "":
        unit = ""
        for i in str(bh_mass):
            if not i.isalnum():
                print(f"{Fore.LIGHTCYAN_EX} mass converter says:{Fore.LIGHTRED_EX} Incorrect mass! cannot have non-alnums!")
                print()
                return False, None
            if i.isalpha():
                unit += str(i)
        if unit != "" and unit != "kg":
            print(f"{Fore.LIGHTCYAN_EX} mass converter says:{Fore.LIGHTRED_EX} Incorrect mass unit! must have 'kg' or without unit!")
            print()
            return False, None
        elif unit != "" and unit == "kg":
            bh_mass = bh_mass.replace("kg", "")
            bh_mass = int(bh_mass)
            solar_mass = bh_mass / sun
        elif unit == "":
            bh_mass = int(bh_mass)
            solar_mass = bh_mass
            bh_mass = bh_mass * sun
            if bh_mass == 0:
                print(f"{Fore.LIGHTCYAN_EX} mass converter says:{Fore.LIGHTRED_EX} mass cannot be zero!")
                print()
                return False, None
        return solar_mass, bh_mass
    else:
        print(f"{Fore.LIGHTCYAN_EX} mass converter says:{Fore.LIGHTRED_EX} mass cannot be empty!")
        print()
        return False, None
    
#selecting which part of the black hole to show
def select_bh_info(bh_type):
    print(Fore.LIGHTYELLOW_EX + "Now Choose which informations of this Black Hole you want to know!\n"
    + Fore.LIGHTBLACK_EX + "Enter 1 to show and 0 to not show." + Fore.RESET)
    print(Fore.LIGHTYELLOW_EX + "==================================================" + Fore.RESET)
    print()
    match bh_type:
        case "schw" | "1":
            #singularity
            while True:
                singularity = input(Fore.LIGHTYELLOW_EX + "Point Singularity?\n"
                    + Fore.LIGHTBLACK_EX + "a point at the very center of the Black Hole which a function takes an infinite value, especially in space–time when matter is infinitely dense!\n"
                    + Fore.RESET + "==============================\n"
                    ">>>").strip().lower()
                if not singularity in ["0", "1"]:
                    print(f"{Fore.LIGHTCYAN_EX} select_bh_info says:{Fore.LIGHTRED_EX} Incorrect entry! must be '1' or '0' !")
                    print()
                    continue
                else:
                    print(f"{Fore.LIGHTGREEN_EX} Sure! {Fore.RESET}")
                    print()
                    if singularity == "1":
                        singularity = True
                    else:
                        singularity = False
                    break

            #inner ergosphere
            inner_ergosphere = False

            #inner event horizon
            inner_event_horizon = False

            #outer event horizon
            while True:
                outer_event_horizon = input(Fore.LIGHTYELLOW_EX + "Outer Event Horizon?\n"
                    + Fore.LIGHTBLACK_EX + "a notional boundary around a black hole beyond which no light or other radiation can escape!\n"
                    + Fore.RESET + "==============================\n"
                    ">>>").strip().lower()
                if not outer_event_horizon in ["0", "1"]:
                    print(f"{Fore.LIGHTCYAN_EX} select_bh_info says:{Fore.LIGHTRED_EX} Incorrect entry! must be '1' or '0' !")
                    print()
                    continue
                else:
                    print(f"{Fore.LIGHTGREEN_EX} Sure! {Fore.RESET}")
                    if outer_event_horizon == "1":
                        outer_event_horizon = True
                    else:
                        outer_event_horizon = False
                    break

            #outer ergosphere
            while True:
                outer_ergosphere = input(Fore.LIGHTYELLOW_EX + "Outer Ergosphere?\n"
                    + Fore.LIGHTBLACK_EX + "the region in between the event horizon and the stationary limit.!\n"
                    + Fore.RESET + "==============================\n"
                    ">>>").strip().lower()
                if not outer_ergosphere in ["0", "1"]:
                    print(f"{Fore.LIGHTCYAN_EX} select_bh_info says:{Fore.LIGHTRED_EX} Incorrect entry! must be '1' or '0' !")
                    print()
                    continue
                else:
                    print(f"{Fore.LIGHTGREEN_EX} Sure! {Fore.RESET}")
                    if outer_ergosphere == "1":
                        outer_ergosphere = True
                    else:
                        outer_ergosphere = False
                    break

            #returning results for schw black hole
            return singularity, inner_ergosphere, inner_event_horizon, outer_event_horizon, outer_ergosphere
        
        case "kerr" | "2":
            #singularity
            while True:
                singularity = input(Fore.LIGHTYELLOW_EX + "Ring Singularity?\n"
                    + Fore.LIGHTBLACK_EX + "A ring singularity or ringularity is the gravitational singularity of a rotating black hole, or a Kerr black hole, that is shaped like a ring which a function takes an infinite value, especially in space–time when matter is infinitely dense!\n"
                    + Fore.RESET + "==============================\n"
                    ">>>").strip().lower()
                if not singularity in ["0", "1"]:
                    print(f"{Fore.LIGHTCYAN_EX} select_bh_info says:{Fore.LIGHTRED_EX} Incorrect entry! must be '1' or '0' !")
                    print()
                    continue
                else:
                    print(f"{Fore.LIGHTGREEN_EX} Sure! {Fore.RESET}")
                    print()
                    if singularity == "1":
                        singularity = True
                    else:
                        singularity = False
                    break

            #inner ergosphere
            while True:
                inner_ergosphere = input(Fore.LIGHTYELLOW_EX + "Inner Ergosphere?\n"
                    + Fore.LIGHTBLACK_EX + "The region between the surface of infinite redshift (outer) and the event horizon (inner)!\n"
                    + Fore.RESET + "==============================\n"
                    ">>>").strip().lower()
                if not inner_ergosphere in ["0", "1"]:
                    print(f"{Fore.LIGHTCYAN_EX} select_bh_info says:{Fore.LIGHTRED_EX} Incorrect entry! must be '1' or '0' !")
                    print()
                    continue
                else:
                    print(f"{Fore.LIGHTGREEN_EX} Sure! {Fore.RESET}")
                    print()
                    if inner_ergosphere == "1":
                        inner_ergosphere = True
                    else:
                        inner_ergosphere = False
                    break

            #inner event horizon
            while True:
                inner_event_horizon = input(Fore.LIGHTYELLOW_EX + "Inner Event Horizon?\n"
                    + Fore.LIGHTBLACK_EX + "The inner ergosphere is the region between the inner event horizon and the inner stationary-limit surface of a rotating Kerr black hole. Inside this region, spacetime is strongly affected by frame dragging, making stationary motion relative to distant observers impossible!\n"
                    + Fore.RESET + "==============================\n"
                    ">>>").strip().lower()
                if not inner_event_horizon in ["0", "1"]:
                    print(f"{Fore.LIGHTCYAN_EX} select_bh_info says:{Fore.LIGHTRED_EX} Incorrect entry! must be '1' or '0' !")
                    print()
                    continue
                else:
                    print(f"{Fore.LIGHTGREEN_EX} Sure! {Fore.RESET}")
                    print()
                    if inner_event_horizon == "1":
                        inner_event_horizon = True
                    else:
                        inner_event_horizon = False
                    break

            #outer event horizon
            while True:
                outer_event_horizon = input(Fore.LIGHTYELLOW_EX + "Outer Event Horizon?\n"
                    + Fore.LIGHTBLACK_EX + "The outer event horizon is the boundary of a black hole beyond which nothing, not even light, can escape to the outside universe. For a rotating Kerr black hole, it is the outermost of the two event horizons and marks the point of no return for infalling objects!\n"
                    + Fore.RESET + "==============================\n"
                    ">>>").strip().lower()
                if not outer_event_horizon in ["0", "1"]:
                    print(f"{Fore.LIGHTCYAN_EX} select_bh_info says:{Fore.LIGHTRED_EX} Incorrect entry! must be '1' or '0' !")
                    print()
                    continue
                else:
                    print(f"{Fore.LIGHTGREEN_EX} Sure! {Fore.RESET}")
                    print()
                    if outer_event_horizon == "1":
                        outer_event_horizon = True
                    else:
                        outer_event_horizon = False
                    break

            #outer ergosphere
            while True:
                outer_ergosphere = input(Fore.LIGHTYELLOW_EX + "Outer Ergosphere?\n"
                    + Fore.LIGHTBLACK_EX + "The outer ergosphere is the region between the outer event horizon and the outer stationary-limit surface of a rotating Kerr black hole. Within this region, frame dragging is so strong that no observer can remain stationary relative to distant space!\n"
                    + Fore.RESET + "==============================\n"
                    ">>>").strip().lower()
                if not outer_ergosphere in ["0", "1"]:
                    print(f"{Fore.LIGHTCYAN_EX} select_bh_info says:{Fore.LIGHTRED_EX} Incorrect entry! must be '1' or '0' !")
                    print()
                    continue
                else:
                    print(f"{Fore.LIGHTGREEN_EX} Sure! {Fore.RESET}")
                    print()
                    if outer_ergosphere == "1":
                        outer_ergosphere = True
                    else:
                        outer_ergosphere = False
                    break

            #returning results for schw black hole
            return singularity, inner_ergosphere, inner_event_horizon, outer_event_horizon, outer_ergosphere
        
def black_hole_info_print(bh_type, bh_mass, solar_mass, a_star, singularity, inner_ergosphere, inner_horizon, outer_horizon, outer_ergosphere):

    print()
    print(Fore.BLUE + "==================================================" + Fore.RESET)
    print(Fore.RED + "Your Black Hole Statistics" + Fore.RESET)

    match bh_type:
        case "schw" | "1":
            #type
            print(Fore.LIGHTGREEN_EX + "Type: " + Fore.RESET + "Non-Spinning Black Hole", end="\n")

            #dimensionless spin parameter
            print(Fore.LIGHTGREEN_EX + f"Dimensionless Spin Parameter: {Fore.RESET} 0", end="\n")

            #printing in solar mass
            print(Fore.LIGHTGREEN_EX + f"Sun Mass: {Fore.RESET} {solar_mass:,} Solar Masses!", end="\n")
            if solar_mass < 1:
                print(Fore.LIGHTBLACK_EX + f"Less than one Solar Mass!" + Fore.RESET)
            else:
                print(Fore.LIGHTBLACK_EX + f"{NumberLex.mag_caller(str(solar_mass))} Solar Masses!" + Fore.RESET)

            #printing in kilograms
            print(Fore.LIGHTGREEN_EX + f"SI Mass{Fore.RESET} {bh_mass:,} KG!", end="\n")
            if bh_mass < 1:
                print(Fore.LIGHTBLACK_EX + f"Less than one KiloGram!" + Fore.RESET)
            else:
                print(Fore.LIGHTBLACK_EX + f"{NumberLex.mag_caller(str(bh_mass))} KiloGrams!" + Fore.RESET)
            print(Fore.RED + "----------------------------------------" + Fore.RESET)

            #singularity
            if singularity == True:
                singularity_radius(bh_mass, a_star, True)
                
            #event horizon
            if outer_horizon == True:
                outer_event_horizon_radius(bh_type, bh_mass, 0, np.radians(0), True)

            #ergosphere
            if outer_ergosphere == True:
                print(Fore.LIGHTGREEN_EX + f"Ergosphere at the Equator: {Fore.RESET} Inside the Horizon!", end="\n")
                print(Fore.RED + "----------------------------------------" + Fore.RESET)

        case "kerr" | "2":
            #type
            print(Fore.LIGHTGREEN_EX + "Type: " + Fore.RESET + "Spinning Black Hole", end="\n")

            #dimensionless spin parameter
            print(Fore.LIGHTGREEN_EX + f"Dimensionless Spin Parameter: {Fore.RESET} {a_star}", end="\n")
            velocity = spin_parser(str(a_star))[1]
            print(Fore.LIGHTGREEN_EX + f"Black Hole Velocity at the Equator: {Fore.RESET} {speed_unit_manager(velocity)}", end="\n")

            #printing in solar mass
            print(Fore.LIGHTGREEN_EX + f"Sun Mass: {Fore.RESET} {solar_mass:,} Solar Masses!", end="\n")
            if solar_mass < 1:
                print(Fore.LIGHTBLACK_EX + f"Less than one Solar Mass!" + Fore.RESET)
            else:
                print(Fore.LIGHTBLACK_EX + f"{NumberLex.mag_caller(str(solar_mass))} Solar Masses!" + Fore.RESET)

            #printing in kilograms
            print(Fore.LIGHTGREEN_EX + f"SI Mass: {Fore.RESET} {bh_mass:,} KG!", end="\n")
            if bh_mass < 1:
                print(Fore.LIGHTBLACK_EX + f"Less than one KiloGram!" + Fore.RESET)
            else:
                print(Fore.LIGHTBLACK_EX + f"{NumberLex.mag_caller(str(bh_mass))} KiloGrams!" + Fore.RESET)
            print(Fore.RED + "----------------------------------------" + Fore.RESET)

            #singularity
            if singularity == True:
                singularity_radius(bh_mass, a_star, True)

            #inner ergosphere
            if inner_ergosphere == True:
                #polar radius
                print(Fore.LIGHTMAGENTA_EX + f"Inner Ergosphere At The Poles:", end="\n")
                inner_ergosphere_radius(bh_mass, a_star, np.radians(0), True)

                #equatorial radius
                print(Fore.LIGHTMAGENTA_EX + f"Inner Ergosphere At The Equator:", end="\n")
                inner_ergosphere_radius(bh_mass, a_star, np.radians(90), True)

            #inner horizon
            if inner_horizon == True:
                #polar radius
                print(Fore.LIGHTMAGENTA_EX + f"Inner Event Horizon At The Poles:", end="\n")
                inner_event_horizon_radius(bh_type, bh_mass, a_star, np.radians(0), True)

                #equatorial radius
                print(Fore.LIGHTMAGENTA_EX + f"Inner Event Horizon At The Equator:", end="\n")
                inner_event_horizon_radius(bh_type, bh_mass, a_star, np.radians(90), True)

            #outer horizon
            if outer_horizon == True:
                #polar radius
                print(Fore.LIGHTMAGENTA_EX + f"Outer Event Horizon At The Poles:", end="\n")
                outer_event_horizon_radius(bh_type, bh_mass, a_star, np.radians(0), True)

                #equatorial radius
                print(Fore.LIGHTMAGENTA_EX + f"Outer Event Horizon At The Equator:", end="\n")
                outer_event_horizon_radius(bh_type, bh_mass, a_star, np.radians(90), True)
                
            #outer ergosphere
            if outer_ergosphere == True:
                #polar radius
                print(Fore.LIGHTMAGENTA_EX + f"Outer Ergosphere At The Poles:", end="\n")
                outer_ergosphere_radius(bh_mass, a_star, np.radians(0), True)

                #equatorial radius
                print(Fore.LIGHTMAGENTA_EX + f"Outer Ergosphere At The Equator:", end="\n")
                outer_ergosphere_radius(bh_mass, a_star, np.radians(90), True)

    print(Fore.BLUE + "==================================================" + Fore.RESET)
    print()

def black_hole_construct(bh_type, bh_mass, a_star, singularity, inner_ergosphere, inner_horizon, outer_horizon, outer_ergosphere):
    match bh_type:
        case "schw" | "1":
            #singularity
            if singularity == True:
                singularity_draw(bh_mass, a_star, singularity_radius(bh_mass, a_star, False), outer_event_horizon_radius(bh_type, bh_mass, a_star, np.radians(90), False))
                
            #event horizon
            if outer_horizon == True:
                outer_event_horizon_draw(bh_mass, bh_type, a_star, outer_event_horizon_radius(bh_type, bh_mass, a_star, np.radians(90), False))

        case "kerr" | "2":
            #singularity
            if singularity == True:
                singularity_draw(bh_mass, a_star, singularity_radius(bh_mass, a_star, False), outer_event_horizon_radius(bh_type, bh_mass, a_star, np.radians(90), False))

            #inner ergosphere
            if inner_ergosphere == True:
                inner_ergosphere_draw(bh_mass, a_star)

            #inner horizon
            if inner_horizon == True:
                inner_event_horizon_draw(bh_type, bh_mass, a_star, inner_event_horizon_radius(bh_type, bh_mass, a_star, np.radians(90), False))

            #inner horizon
            if outer_horizon == True:
                outer_event_horizon_draw(bh_mass, bh_type, a_star, outer_event_horizon_radius(bh_type, bh_mass, a_star, np.radians(90), False))

            #outer ergosphere
            if outer_ergosphere == True:
                outer_ergosphere_draw(bh_type, bh_mass, a_star)

#selecting the type of space-time to show
def select_space_time_type(bh_type):
    exit_first_loop = 0
    dimension, flat_st, bent_st, twisted_st = False, False, False, False

    print(Fore.LIGHTYELLOW_EX + "Now Choose which type of space-time you want to show!\n"
    + Fore.LIGHTBLACK_EX + "Enter 1 to show and 0 to not show." + Fore.RESET)
    print(Fore.LIGHTYELLOW_EX + "==================================================" + Fore.RESET)
    print()

    while True:
        dimension = input(Fore.LIGHTYELLOW_EX + "2D Space-Time or 3D Space_Time?\n"
            + Fore.LIGHTBLACK_EX + "2D Space-Time: A simplified visualization using one spatial "
            "dimension and one time dimension, allowing the curvature of "
            "spacetime to be represented on a two-dimensional diagram.\n"
            "3D Space-Time: An extended visualization using two spatial "
            "dimensions and one time dimension, providing a more detailed "
            "representation of the geometry and curvature of spacetime.\n"
            "Enter 1 to show 2D and 2 to show 3D or 0 to not show.\n"
            + Fore.RESET + "==============================\n"
            ">>>").strip().lower()

        if not dimension in ["0", "1", "2"]:
            print(f"{Fore.LIGHTCYAN_EX} select_space_time_type says:{Fore.LIGHTRED_EX} Incorrect entry! must be '1' or '2' !")
            print()
            continue
        else:
            print(f"{Fore.LIGHTGREEN_EX} Sure! {Fore.RESET}")
            print()
            match dimension:
                case "0":
                    return False, False, False, False
                case "1":
                    dimension = "2D"
                case "2":
                    dimension = "3D"
            break

    #flat space-time
    while True:
        flat_st = input(Fore.LIGHTYELLOW_EX + "Flat Space-time?\n"
            + Fore.LIGHTBLACK_EX + "A flat spacetime is a region of spacetime where there is no significant curvature caused by gravity. It is described by Minkowski spacetime, in which space and time form a four-dimensional geometry with a constant metric.\n"
            + Fore.RESET + "==============================\n"
            ">>>").strip().lower()
        if not flat_st in ["0", "1"]:
            print(f"{Fore.LIGHTCYAN_EX} select_space_time_type says:{Fore.LIGHTRED_EX} Incorrect entry! must be '1' or '0' !")
            print()
            continue
        else:
            print(f"{Fore.LIGHTGREEN_EX} Sure! {Fore.RESET}")
            print()
            if flat_st == "1":
                ax.grid(False)
                flat_st = True
                #bent space-time
                while True:
                    bent_st = input(Fore.LIGHTYELLOW_EX + "Bend the Space-time?\n"
                    + Fore.LIGHTBLACK_EX + "A curved spacetime is a region of spacetime where matter and energy distort the geometry of space and time. It is described by general relativity, in which gravity is interpreted as the curvature of spacetime rather than a conventional force.\n"
                    + Fore.RESET + "==============================\n"
                    ">>>").strip().lower()
                    if not bent_st in ["0", "1"]:
                        print(f"{Fore.LIGHTCYAN_EX} select_space_time_type says:{Fore.LIGHTRED_EX} Incorrect entry! must be '1' or '0' !")
                        print()
                        continue
                    else:
                        print(f"{Fore.LIGHTGREEN_EX} Sure! {Fore.RESET}")
                        print()
                        if bent_st == "1":
                            bent_st = True
                        else:
                            bent_st = False

                        exit_first_loop = 1
                        break
                        
                #twisted space-time
                if bh_type == "kerr" or bh_type == "2":
                    while True:
                        twisted_st = input(Fore.LIGHTYELLOW_EX + "Twist the Space-time?\n"
                        + Fore.LIGHTBLACK_EX + "A twisted spacetime refers to spacetime whose geometry is affected not only by curvature, but also by the rotation of a massive object. In general relativity, rotating matter can drag nearby spacetime around with it, producing frame dragging!\n"
                        + Fore.RESET + "==============================\n"
                        ">>>").strip().lower()
                        if not twisted_st in ["0", "1"]:
                            print(f"{Fore.LIGHTCYAN_EX} select_space_time_type says:{Fore.LIGHTRED_EX} Incorrect entry! must be '1' or '0' !")
                            print()
                            continue
                        else:
                            print(f"{Fore.LIGHTGREEN_EX} Sure! {Fore.RESET}")
                            print()
                            if twisted_st == "1":
                                twisted_st = True
                            else:
                                twisted_st = False

                            exit_first_loop = 1
                            break

                if exit_first_loop == 1:
                    break
            else:
                flat_st = False
                break

    return dimension, flat_st, bent_st, twisted_st

def space_time_construct(bh_type, bh_mass, a_star, dimension, flat_st, bent_st, twisted_st, event_horizon_radius):
    #flat space-time
    if flat_st == True:
        space_time = space_time_draw(dimension, event_horizon_radius)

    #bent space-time
    if bent_st == True:
        bent_space_time = bend_space_time(space_time[0], space_time[1], space_time[2], space_time[3], space_time[4], dimension, event_horizon_radius)

    #twisted space-time
    if (bh_type == "kerr" or bh_type == "2") and (twisted_st == True):
        if bent_st == True:
            twisted_space_time = twist_space_time(bh_mass, a_star, bent_space_time[0], bent_space_time[1], bent_space_time[2], bent_space_time[3], bent_space_time[4], dimension, event_horizon_radius)
        else:
            twisted_space_time = twist_space_time(bh_mass, a_star, space_time[0], space_time[1], space_time[2], space_time[3], space_time[4], dimension, event_horizon_radius)

#selecting photon orbitals to show
def select_photon_orbitals(bh_type):
    exit_first_loop = 0
    co_photon_ring, co_photon_multiple_rings, co_photon_sphere, counter_photon_ring, counter_photon_multiple_rings, counter_photon_sphere = False, False, False, False, False, False
    co_photon_ring_angle, co_photon_multiple_rings_count, counter_photon_ring_angle, counter_photon_multiple_rings_count = 0, 0, 0, 0
    
    print(Fore.LIGHTYELLOW_EX + "Now Choose which Photon orbitals of this Black Hole you want to show!\n"
    + Fore.LIGHTBLACK_EX + "Enter 1 to show and 0 to not show." + Fore.RESET)
    print(Fore.LIGHTYELLOW_EX + "==================================================" + Fore.RESET)
    print()

    match bh_type:
        case "schw" | "1":
            #photon sphere
            while True:
                photon_sphere = input(Fore.LIGHTYELLOW_EX + "Photon Sphere?\n"
                + Fore.LIGHTBLACK_EX + "A photon sphere is a spherical region of space surrounding a very compact object — like a black hole — where gravity is so intense that light (photons) is forced to travel in a circular orbit!\n"
                "Enter 1 to show and 0 to not show.\n"
                + Fore.RESET + "==============================\n"
                ">>>").strip().lower()

                if not photon_sphere in ["0", "1"]:
                    print(f"{Fore.LIGHTCYAN_EX} select_photon_orbital says:{Fore.LIGHTRED_EX} Incorrect entry! must be '1' or '0'!")
                    print()
                    continue
                if photon_sphere == "0": break
                while True:
                    photon_sphere = input(Fore.LIGHTYELLOW_EX + "Choose Type?\n"
                    + Fore.LIGHTBLACK_EX + 
                    "1 for a Single Photon Ring\n"
                    "2 for Multiple Rings\n"
                    "3 for Full Spherical\n"
                    + Fore.RESET + "==============================\n"
                    ">>>").strip().lower()
                    if not photon_sphere in ["1", "2", "3"]:
                        print(f"{Fore.LIGHTCYAN_EX} select_photon_orbital says:{Fore.LIGHTRED_EX} Incorrect Entry! must be '1' or '2' or '3'!")
                        print()
                        continue
                    match photon_sphere:
                        case "1":
                            co_photon_ring = True
                            while True:
                                try:
                                    co_photon_ring_angle = int(input(Fore.LIGHTYELLOW_EX + "Enter the Angle.\n"
                                    + Fore.LIGHTBLACK_EX + "must be between 0° and 180°.\n"
                                    + Fore.RESET + "==============================\n"
                                    ">>>"))
                                    if not 0 <= co_photon_ring_angle <= 180:
                                        print(f"{Fore.LIGHTCYAN_EX} select_photon_orbital says:{Fore.LIGHTRED_EX} Incorrect Angle! must be between '0' and '180'!")
                                        print()
                                        continue
                                    else:
                                        print(f"{Fore.LIGHTGREEN_EX} Photon Ring will draw at the Angle of {co_photon_ring_angle}°! {Fore.RESET}")
                                        print()
                                    break
                                except ValueError:
                                        print(f"{Fore.LIGHTCYAN_EX} select_photon_orbital says:{Fore.LIGHTRED_EX} Incorrect Angle! must be between '0' and '180'!")
                                        print()
                                        continue

                        case "2":
                            co_photon_multiple_rings = True
                            while True:
                                try:
                                    co_photon_multiple_rings_count = int(input(Fore.LIGHTYELLOW_EX + "How Many Rings?\n"
                                    + Fore.LIGHTBLACK_EX + "Enter count between 1(for Equator) to 9.\n"
                                    + Fore.RESET + "==============================\n"
                                    ">>>"))
                                    if not 1 <= co_photon_multiple_rings_count <= 9:
                                        print(f"{Fore.LIGHTCYAN_EX} select_photon_orbital says:{Fore.LIGHTRED_EX} Incorrect Rings Count! must be between '1' and '9'!")
                                        print()
                                        continue
                                    else:
                                        print(f"{Fore.LIGHTGREEN_EX} {co_photon_multiple_rings_count} Photon Rings will be shown! {Fore.RESET}")
                                        print()
                                    break
                                except ValueError:
                                        print(f"{Fore.LIGHTCYAN_EX} select_photon_orbital says:{Fore.LIGHTRED_EX} Incorrect Rings Count! must be between '1' and '9'!")
                                        print()
                                        continue

                        case "3":
                            co_photon_sphere = True
                            print(f"{Fore.LIGHTGREEN_EX} Photon Sphere will be shown as Full Spherical! {Fore.RESET}")
                            print()
                    exit_first_loop = 1
                    break
                if exit_first_loop == 1:
                    break

        case "kerr" | "2":
            #co-rotating photon ring
            while True:
                photon_sphere = input(Fore.LIGHTYELLOW_EX + "Co-Rotating Photon Sphere?\n"
                + Fore.LIGHTBLACK_EX + "A photon orbit that moves in the same direction as the black hole's rotation. Due to frame dragging, the black hole's spin pulls spacetime around with it, allowing photons to orbit closer to the event horizon.\n"
                "Enter 1 to show and 0 to not show.\n"
                + Fore.RESET + "==============================\n"
                ">>>").strip().lower()
                if not photon_sphere in ["0", "1"]:
                    print(f"{Fore.LIGHTCYAN_EX} select_photon_orbital says:{Fore.LIGHTRED_EX} Incorrect entry! must be '1' or '0' !")
                    print()
                    continue
                if photon_sphere == "0": break
                while True:
                    photon_sphere = input(Fore.LIGHTYELLOW_EX + "Choose Type?\n"
                    + Fore.LIGHTBLACK_EX + 
                    "1 for a Single Co-Photon Ring\n"
                    "2 for Multiple Co-Photon Rings\n"
                    "3 for Full Spherical Co-Photon\n"
                    + Fore.RESET + "==============================\n"
                    ">>>").strip().lower()
                    if not photon_sphere in ["1", "2", "3"]:
                        print(f"{Fore.LIGHTCYAN_EX} select_photon_orbital says:{Fore.LIGHTRED_EX} Incorrect Entry! must be '1' or '2' or '3'!")
                        print()
                        continue
                    match photon_sphere:
                        case "1":
                            co_photon_ring = True
                            while True:
                                try:
                                    co_photon_ring_angle = int(input(Fore.LIGHTYELLOW_EX + "Enter the Co-Photon Angle.\n"
                                    + Fore.LIGHTBLACK_EX + "must be between 0° and 180°.\n"
                                    + Fore.RESET + "==============================\n"
                                    ">>>"))
                                    if not 0 <= co_photon_ring_angle <= 180:
                                        print(f"{Fore.LIGHTCYAN_EX} select_photon_orbital says:{Fore.LIGHTRED_EX} Incorrect Co-Photon Angle! must be between '0' and '180'!")
                                        print()
                                        continue
                                    else:
                                        print(f"{Fore.LIGHTGREEN_EX} Co-Photon Ring will draw at the Angle of {co_photon_ring_angle}°! {Fore.RESET}")
                                        print()
                                    break
                                except ValueError:
                                        print(f"{Fore.LIGHTCYAN_EX} select_photon_orbital says:{Fore.LIGHTRED_EX} Incorrect Co-Photon Angle! must be between '0' and '180'!")
                                        print()
                                        continue

                        case "2":
                            co_photon_multiple_rings = True
                            while True:
                                try:
                                    co_photon_multiple_rings_count = int(input(Fore.LIGHTYELLOW_EX + "How Many Co-Photon Rings?\n"
                                    + Fore.LIGHTBLACK_EX + "Enter count between 1(for Equator) to 9.\n"
                                    + Fore.RESET + "==============================\n"
                                    ">>>"))
                                    if not 1 <= co_photon_multiple_rings_count <= 9:
                                        print(f"{Fore.LIGHTCYAN_EX} select_photon_orbital says:{Fore.LIGHTRED_EX} Incorrect Co-Photon Rings Count! must be between '1' and '9'!")
                                        print()
                                        continue
                                    else:
                                        print(f"{Fore.LIGHTGREEN_EX} {co_photon_multiple_rings_count} Co-Photon Rings will be shown! {Fore.RESET}")
                                        print()
                                    break
                                except ValueError:
                                        print(f"{Fore.LIGHTCYAN_EX} select_photon_orbital says:{Fore.LIGHTRED_EX} Incorrect Co-Photon Rings Count! must be between '1' and '9'!")
                                        print()
                                        continue

                        case "3":
                            co_photon_sphere = True
                            print(f"{Fore.LIGHTGREEN_EX} Co-Photon Sphere will be shown as Full Spherical! {Fore.RESET}")
                            print()
                    exit_first_loop = 1
                    break
                if exit_first_loop == 1:
                    break                    
            #counter-rotating photon ring
            while True:
                photon_sphere = input(Fore.LIGHTYELLOW_EX + "Counter-Rotating Photon Sphere?\n"
                + Fore.LIGHTBLACK_EX + "The counter-rotating photon sphere is the unstable circular orbit of light around a rotating black hole where the photon travels opposite to the black hole's spin. Due to frame dragging, it lies farther from the black hole than the co-rotating photon sphere.\n"
                "Enter 1 to show and 0 to not show.\n"
                + Fore.RESET + "==============================\n"
                ">>>").strip().lower()
                if not photon_sphere in ["0", "1"]:
                    print(f"{Fore.LIGHTCYAN_EX} select_photon_orbital says:{Fore.LIGHTRED_EX} Incorrect entry! must be '1' or '0' !")
                    print()
                    continue
                if photon_sphere == "0": break
                while True:
                    photon_sphere = input(Fore.LIGHTYELLOW_EX + "Choose Type?\n"
                    + Fore.LIGHTBLACK_EX + 
                    "1 for a Single Counter-Photon Ring\n"
                    "2 for Multiple Counter-Photon Rings\n"
                    "3 for Full Spherical Counter-Photon\n"
                    + Fore.RESET + "==============================\n"
                    ">>>").strip().lower()
                    if not photon_sphere in ["1", "2", "3"]:
                        print(f"{Fore.LIGHTCYAN_EX} select_photon_orbital says:{Fore.LIGHTRED_EX} Incorrect Entry! must be '1' or '2' or '3'!")
                        print()
                        continue
                    match photon_sphere:
                        case "1":
                            counter_photon_ring = True
                            while True:
                                try:
                                    counter_photon_ring_angle = int(input(Fore.LIGHTYELLOW_EX + "Enter the Counter-Photon Angle.\n"
                                    + Fore.LIGHTBLACK_EX + "must be between 0° and 180°.\n"
                                    + Fore.RESET + "==============================\n"
                                    ">>>"))
                                    if not 0 <= counter_photon_ring_angle <= 180:
                                        print(f"{Fore.LIGHTCYAN_EX} select_photon_orbital says:{Fore.LIGHTRED_EX} Incorrect Counter-Photon Angle! must be between '0' and '180'!")
                                        print()
                                        continue
                                    else:
                                        print(f"{Fore.LIGHTGREEN_EX} Counter-Photon Ring will draw at the Angle of {counter_photon_ring_angle}°! {Fore.RESET}")
                                        print()
                                    break
                                except ValueError:
                                        print(f"{Fore.LIGHTCYAN_EX} select_photon_orbital says:{Fore.LIGHTRED_EX} Incorrect Counter-Photon Angle! must be between '0' and '180'!")
                                        print()
                                        continue

                        case "2":
                            counter_photon_multiple_rings = True
                            while True:
                                try:
                                    counter_photon_multiple_rings_count = int(input(Fore.LIGHTYELLOW_EX + "How Many Counter-Photon Rings?\n"
                                    + Fore.LIGHTBLACK_EX + "Enter count between 1(for Equator) to 9.\n"
                                    + Fore.RESET + "==============================\n"
                                    ">>>"))
                                    if not 1 <= counter_photon_multiple_rings_count <= 9:
                                        print(f"{Fore.LIGHTCYAN_EX} select_photon_orbital says:{Fore.LIGHTRED_EX} Incorrect Counter-Photon Rings Count! must be between '1' and '9'!")
                                        print()
                                        continue
                                    else:
                                        print(f"{Fore.LIGHTGREEN_EX} {counter_photon_multiple_rings_count} Counter-Photon Rings will be shown! {Fore.RESET}")
                                        print()
                                    break
                                except ValueError:
                                    print(f"{Fore.LIGHTCYAN_EX} select_photon_orbital says:{Fore.LIGHTRED_EX} Incorrect Counter-Photon Rings Count! must be between '1' and '9'!")
                                    print()
                                    continue

                        case "3":
                            counter_photon_sphere = True
                            print(f"{Fore.LIGHTGREEN_EX} Counter-Photon Sphere will be shown as Full Spherical! {Fore.RESET}")
                            print()
                    exit_first_loop = 1
                    break
                if exit_first_loop == 1:
                    break      
    return (
        co_photon_ring, 
        co_photon_ring_angle, 
        co_photon_multiple_rings, 
        co_photon_multiple_rings_count, 
        co_photon_sphere,
        counter_photon_ring, 
        counter_photon_ring_angle, 
        counter_photon_multiple_rings, 
        counter_photon_multiple_rings_count, 
        counter_photon_sphere
    )

def photon_orbitals_info_print(bh_type, bh_mass, a_star, co_photon_ring, co_photon_ring_angle, co_photon_multiple_rings, co_photon_multiple_rings_count, co_photon_sphere, counter_photon_ring, counter_photon_ring_angle, counter_photon_multiple_rings, counter_photon_multiple_rings_count, counter_photon_sphere):
    print()
    print(Fore.BLUE + "==================================================" + Fore.RESET)
    print(Fore.RED + "Your Black Hole Photon Orbitals" + Fore.RESET)

    match bh_type:
        case "schw" | "1":
            #single photon ring
            if co_photon_ring == True:
                print(Fore.LIGHTMAGENTA_EX + f"Single Photon Ring:", end="\n")

                photon_sphere_radius(bh_type, bh_mass, a_star, "", np.radians(co_photon_ring_angle), True)

            #multiple photon rings
            elif co_photon_multiple_rings == True:
                print(Fore.LIGHTMAGENTA_EX + f"Multiple Photon Rings:", end="\n")

                additional_angle = 180 / co_photon_multiple_rings_count
                angle = 0
                for i in range(co_photon_multiple_rings_count):
                    photon_sphere_radius(bh_type, bh_mass, a_star, "", np.radians(angle), True)
                    angle += additional_angle

            #photon sphere
            elif co_photon_sphere == True:
                print(Fore.LIGHTMAGENTA_EX + f"Photon Sphere:", end="\n")

                additional_angle = 180 / 9
                angle = 0
                for i in range(9):
                    photon_sphere_radius(bh_type, bh_mass, a_star, "", np.radians(angle), True)
                    angle += additional_angle

        case "kerr" | "2":
            #co-rotating single photon ring
            if co_photon_ring == True:
                print(Fore.LIGHTMAGENTA_EX + f"Co-Rotating Single Photon Ring:", end="\n")

                photon_sphere_radius(bh_type, bh_mass, a_star, "prograde", np.radians(co_photon_ring_angle), True)

            #co-rotating multiple photon rings
            elif co_photon_multiple_rings == True:
                print(Fore.LIGHTMAGENTA_EX + f"Co-Rotating Multiple Photon Rings:", end="\n")

                additional_angle = 180 / co_photon_multiple_rings_count
                angle = 0
                for i in range(co_photon_multiple_rings_count):
                    photon_sphere_radius(bh_type, bh_mass, a_star, "prograde", np.radians(angle), True)
                    angle += additional_angle

            #co-rotating photon sphere
            elif co_photon_sphere == True:
                print(Fore.LIGHTMAGENTA_EX + f"Co-Rotating Photon Sphere:", end="\n")

                additional_angle = 180 / 9
                angle = 0
                for i in range(9):
                    photon_sphere_radius(bh_type, bh_mass, a_star, "prograde", np.radians(angle), True)
                    angle += additional_angle

            #counter-rotating single photon ring
            if counter_photon_ring == True:
                print(Fore.LIGHTMAGENTA_EX + f"Counter-Rotating Single Photon Ring:", end="\n")

                photon_sphere_radius(bh_type, bh_mass, a_star, "retrograde", np.radians(counter_photon_ring_angle), True)

            #counter-rotating photon rings
            if counter_photon_multiple_rings == True:

                print(Fore.LIGHTMAGENTA_EX + f"Counter-Rotating Multiple Photon Rings:", end="\n")

                additional_angle = 180 / counter_photon_multiple_rings_count
                angle = 0
                for i in range(counter_photon_multiple_rings_count):
                    photon_sphere_radius(bh_type, bh_mass, a_star, "retrograde", np.radians(angle), True)
                    angle += additional_angle

            #co-rotating photon sphere
            elif counter_photon_sphere == True:
                print(Fore.LIGHTMAGENTA_EX + f"Counter-Rotating Photon Sphere:", end="\n")

                additional_angle = 180 / 9
                angle = 0
                for i in range(9):
                    photon_sphere_radius(bh_type, bh_mass, a_star, "retrograde", np.radians(angle), True)
                    angle += additional_angle

def photon_orbitals_construct(bh_type, bh_mass, a_star, co_photon_ring, co_photon_ring_angle, co_photon_multiple_rings, co_photon_multiple_rings_count, co_photon_sphere, counter_photon_ring, counter_photon_ring_angle, counter_photon_multiple_rings, counter_photon_multiple_rings_count, counter_photon_sphere):
    match bh_type:
        case "schw" | "1":
            #single photon ring
            if co_photon_ring == True:
                photon_ring_draw(bh_type, bh_mass, a_star, "", np.radians(co_photon_ring_angle))
            
            #multiple photon rings
            elif co_photon_multiple_rings == True:
                additional_angle = 180 / co_photon_multiple_rings_count
                angle = 0
                for i in range(co_photon_multiple_rings_count):
                    photon_ring_draw(bh_type, bh_mass, a_star, "", np.radians(angle))
                    angle += additional_angle

            #photon sphere
            elif co_photon_sphere == True:
                photon_sphere_draw(bh_type, bh_mass, a_star, "")
                
        case "kerr" | "2":
            #co-rotating single photon ring
            if co_photon_ring == True:
                photon_ring_draw(bh_type, bh_mass, a_star, "prograde", np.radians(co_photon_ring_angle))

            #co-rotating multiple photon rings
            elif co_photon_multiple_rings == True:
                additional_angle = 180 / co_photon_multiple_rings_count
                angle = 0
                for i in range(co_photon_multiple_rings_count):
                    photon_ring_draw(bh_type, bh_mass, a_star, "prograde", np.radians(angle))
                    angle += additional_angle

            #co-rotating photon sphere
            elif co_photon_sphere == True:
                photon_sphere_draw(bh_type, bh_mass, a_star, "prograde")

            #counter-rotating single photon ring
            if counter_photon_ring == True:
                photon_ring_draw(bh_type, bh_mass, a_star, "retrograde", np.radians(counter_photon_ring_angle))


            #counter-rotating multiple photon rings
            elif counter_photon_multiple_rings == True:
                additional_angle = 180 / counter_photon_multiple_rings_count
                angle = 0
                for i in range(counter_photon_multiple_rings_count):
                    photon_ring_draw(bh_type, bh_mass, a_star, "retrograde", np.radians(angle))
                    angle += additional_angle

            #counter-rotating photon sphere
            elif counter_photon_sphere == True:
                photon_sphere_draw(bh_type, bh_mass, a_star, "retrograde")

#selecting ISCO orbitals to show
def select_isco_orbitals(bh_type):
    exit_first_loop = 0
    co_isco_ring, co_isco_multiple_rings, co_isco_sphere, counter_isco_ring, counter_isco_multiple_rings, counter_isco_sphere = False, False, False, False, False, False
    co_isco_ring_angle, co_isco_multiple_rings_count, counter_isco_ring_angle, counter_isco_multiple_rings_count = 0, 0, 0, 0
    
    print(Fore.LIGHTYELLOW_EX + "Now Choose which ISCO orbitals of this Black Hole you want to show!\n"
    + Fore.LIGHTBLACK_EX + "Enter 1 to show and 0 to not show." + Fore.RESET)
    print(Fore.LIGHTYELLOW_EX + "==================================================" + Fore.RESET)
    print()

    match bh_type:
        case "schw" | "1":
            #isco sphere
            while True:
                isco_sphere = input(Fore.LIGHTYELLOW_EX + "ISCO Sphere?\n"
                + Fore.LIGHTBLACK_EX + "The innermost stable circular orbit (often called the ISCO) is the smallest marginally stable circular orbit in which a test particle can stably orbit a massive object in general relativity.\n"
                "Enter 1 to show and 0 to not show.\n"
                + Fore.RESET + "==============================\n"
                ">>>").strip().lower()

                if not isco_sphere in ["0", "1"]:
                    print(f"{Fore.LIGHTCYAN_EX} select_isco_orbital says:{Fore.LIGHTRED_EX} Incorrect entry! must be '1' or '0'!")
                    print()
                    continue
                if isco_sphere == "0": break
                while True:
                    isco_sphere = input(Fore.LIGHTYELLOW_EX + "Choose Type?\n"
                    + Fore.LIGHTBLACK_EX + 
                    "1 for a Single ISCO Ring\n"
                    "2 for Multiple Rings\n"
                    "3 for Full Spherical\n"
                    + Fore.RESET + "==============================\n"
                    ">>>").strip().lower()
                    if not isco_sphere in ["1", "2", "3"]:
                        print(f"{Fore.LIGHTCYAN_EX} select_isco_orbital says:{Fore.LIGHTRED_EX} Incorrect Entry! must be '1' or '2' or '3'!")
                        print()
                        continue
                    match isco_sphere:
                        case "1":
                            co_isco_ring = True
                            while True:
                                try:
                                    co_isco_ring_angle = int(input(Fore.LIGHTYELLOW_EX + "Enter the Angle.\n"
                                    + Fore.LIGHTBLACK_EX + "must be between 0° and 180°.\n"
                                    + Fore.RESET + "==============================\n"
                                    ">>>"))
                                    if not 0 <= co_isco_ring_angle <= 180:
                                        print(f"{Fore.LIGHTCYAN_EX} select_isco_orbital says:{Fore.LIGHTRED_EX} Incorrect Angle! must be between '0' and '180'!")
                                        print()
                                        continue
                                    else:
                                        print(f"{Fore.LIGHTGREEN_EX} ISCO Ring will draw at the Angle of {co_isco_ring_angle}°! {Fore.RESET}")
                                        print()
                                    break
                                except ValueError:
                                        print(f"{Fore.LIGHTCYAN_EX} select_isco_orbital says:{Fore.LIGHTRED_EX} Incorrect Angle! must be between '0' and '180'!")
                                        print()
                                        continue

                        case "2":
                            co_isco_multiple_rings = True
                            while True:
                                try:
                                    co_isco_multiple_rings_count = int(input(Fore.LIGHTYELLOW_EX + "How Many Rings?\n"
                                    + Fore.LIGHTBLACK_EX + "Enter count between 1(for Equator) to 9.\n"
                                    + Fore.RESET + "==============================\n"
                                    ">>>"))
                                    if not 1 <= co_isco_multiple_rings_count <= 9:
                                        print(f"{Fore.LIGHTCYAN_EX} select_isco_orbital says:{Fore.LIGHTRED_EX} Incorrect Rings Count! must be between '1' and '9'!")
                                        print()
                                        continue
                                    else:
                                        print(f"{Fore.LIGHTGREEN_EX} {co_isco_multiple_rings_count} ISCO Rings will be shown! {Fore.RESET}")
                                        print()
                                    break
                                except ValueError:
                                        print(f"{Fore.LIGHTCYAN_EX} select_isco_orbital says:{Fore.LIGHTRED_EX} Incorrect Rings Count! must be between '1' and '9'!")
                                        print()
                                        continue

                        case "3":
                            co_isco_sphere = True
                            print(f"{Fore.LIGHTGREEN_EX} ISCO Sphere will be shown as Full Spherical! {Fore.RESET}")
                            print()
                    exit_first_loop = 1
                    break
                if exit_first_loop == 1:
                    break

        case "kerr" | "2":
            #co-rotating isco ring
            while True:
                isco_sphere = input(Fore.LIGHTYELLOW_EX + "Co-Rotating ISCO Sphere?\n"
                + Fore.LIGHTBLACK_EX + "an ISCO orbit that moves in the same direction as the black hole's rotation. Due to frame dragging, the black hole's spin pulls spacetime around with it, allowing matter to orbit closer to the event horizon.\n"
                "Enter 1 to show and 0 to not show.\n"
                + Fore.RESET + "==============================\n"
                ">>>").strip().lower()
                if not isco_sphere in ["0", "1"]:
                    print(f"{Fore.LIGHTCYAN_EX} select_isco_orbital says:{Fore.LIGHTRED_EX} Incorrect entry! must be '1' or '0' !")
                    print()
                    continue
                if isco_sphere == "0": break
                while True:
                    isco_sphere = input(Fore.LIGHTYELLOW_EX + "Choose Type?\n"
                    + Fore.LIGHTBLACK_EX + 
                    "1 for a Single Co-ISCO Ring\n"
                    "2 for Multiple Co-ISCO Rings\n"
                    "3 for Full Spherical Co-ISCO\n"
                    + Fore.RESET + "==============================\n"
                    ">>>").strip().lower()
                    if not isco_sphere in ["1", "2", "3"]:
                        print(f"{Fore.LIGHTCYAN_EX} select_isco_orbital says:{Fore.LIGHTRED_EX} Incorrect Entry! must be '1' or '2' or '3'!")
                        print()
                        continue
                    match isco_sphere:
                        case "1":
                            co_isco_ring = True
                            while True:
                                try:
                                    co_isco_ring_angle = int(input(Fore.LIGHTYELLOW_EX + "Enter the Co-ISCO Angle.\n"
                                    + Fore.LIGHTBLACK_EX + "must be between 0° and 180°.\n"
                                    + Fore.RESET + "==============================\n"
                                    ">>>"))
                                    if not 0 <= co_isco_ring_angle <= 180:
                                        print(f"{Fore.LIGHTCYAN_EX} select_isco_orbital says:{Fore.LIGHTRED_EX} Incorrect Co-ISCO Angle! must be between '0' and '180'!")
                                        print()
                                        continue
                                    else:
                                        print(f"{Fore.LIGHTGREEN_EX} Co-ISCO Ring will draw at the Angle of {co_isco_ring_angle}°! {Fore.RESET}")
                                        print()
                                    break
                                except ValueError:
                                        print(f"{Fore.LIGHTCYAN_EX} select_isco_orbital says:{Fore.LIGHTRED_EX} Incorrect Co-ISCO Angle! must be between '0' and '180'!")
                                        print()
                                        continue

                        case "2":
                            co_isco_multiple_rings = True
                            while True:
                                try:
                                    co_isco_multiple_rings_count = int(input(Fore.LIGHTYELLOW_EX + "How Many Co-ISCO Rings?\n"
                                    + Fore.LIGHTBLACK_EX + "Enter count between 1(for Equator) to 9.\n"
                                    + Fore.RESET + "==============================\n"
                                    ">>>"))
                                    if not 1 <= co_isco_multiple_rings_count <= 9:
                                        print(f"{Fore.LIGHTCYAN_EX} select_isco_orbital says:{Fore.LIGHTRED_EX} Incorrect Co-ISCO Rings Count! must be between '1' and '9'!")
                                        print()
                                        continue
                                    else:
                                        print(f"{Fore.LIGHTGREEN_EX} {co_isco_multiple_rings_count} Co-ISCO Rings will be shown! {Fore.RESET}")
                                        print()
                                    break
                                except ValueError:
                                        print(f"{Fore.LIGHTCYAN_EX} select_isco_orbital says:{Fore.LIGHTRED_EX} Incorrect Co-ISCO Rings Count! must be between '1' and '9'!")
                                        print()
                                        continue

                        case "3":
                            co_isco_sphere = True
                            print(f"{Fore.LIGHTGREEN_EX} Co-ISCO Sphere will be shown as Full Spherical! {Fore.RESET}")
                            print()
                    exit_first_loop = 1
                    break
                if exit_first_loop == 1:
                    break                    
            #counter-rotating isco ring
            while True:
                isco_sphere = input(Fore.LIGHTYELLOW_EX + "Counter-Rotating ISCO Sphere?\n"
                + Fore.LIGHTBLACK_EX + "The counter-rotating ISCO sphere is the unstable circular orbit of matter around a rotating black hole where the matter travels opposite to the black hole's spin. Due to frame dragging, it lies farther from the black hole than the co-rotating ISCO sphere.\n"
                "Enter 1 to show and 0 to not show.\n"
                + Fore.RESET + "==============================\n"
                ">>>").strip().lower()
                if not isco_sphere in ["0", "1"]:
                    print(f"{Fore.LIGHTCYAN_EX} select_isco_orbital says:{Fore.LIGHTRED_EX} Incorrect entry! must be '1' or '0' !")
                    print()
                    continue
                if isco_sphere == "0": break
                while True:
                    isco_sphere = input(Fore.LIGHTYELLOW_EX + "Choose Type?\n"
                    + Fore.LIGHTBLACK_EX + 
                    "1 for a Single Counter-ISCO Ring\n"
                    "2 for Multiple Counter-ISCO Rings\n"
                    "3 for Full Spherical Counter-ISCO\n"
                    + Fore.RESET + "==============================\n"
                    ">>>").strip().lower()
                    if not isco_sphere in ["1", "2", "3"]:
                        print(f"{Fore.LIGHTCYAN_EX} select_isco_orbital says:{Fore.LIGHTRED_EX} Incorrect Entry! must be '1' or '2' or '3'!")
                        print()
                        continue
                    match isco_sphere:
                        case "1":
                            counter_isco_ring = True
                            while True:
                                try:
                                    counter_isco_ring_angle = int(input(Fore.LIGHTYELLOW_EX + "Enter the Counter-ISCO Angle.\n"
                                    + Fore.LIGHTBLACK_EX + "must be between 0° and 180°.\n"
                                    + Fore.RESET + "==============================\n"
                                    ">>>"))
                                    if not 0 <= counter_isco_ring_angle <= 180:
                                        print(f"{Fore.LIGHTCYAN_EX} select_isco_orbital says:{Fore.LIGHTRED_EX} Incorrect Counter-ISCO Angle! must be between '0' and '180'!")
                                        print()
                                        continue
                                    else:
                                        print(f"{Fore.LIGHTGREEN_EX} Counter-ISCO Ring will draw at the Angle of {counter_isco_ring_angle}°! {Fore.RESET}")
                                        print()
                                    break
                                except ValueError:
                                        print(f"{Fore.LIGHTCYAN_EX} select_isco_orbital says:{Fore.LIGHTRED_EX} Incorrect Counter-ISCO Angle! must be between '0' and '180'!")
                                        print()
                                        continue

                        case "2":
                            counter_isco_multiple_rings = True
                            while True:
                                try:
                                    counter_isco_multiple_rings_count = int(input(Fore.LIGHTYELLOW_EX + "How Many Counter-ISCO Rings?\n"
                                    + Fore.LIGHTBLACK_EX + "Enter count between 1(for Equator) to 9.\n"
                                    + Fore.RESET + "==============================\n"
                                    ">>>"))
                                    if not 1 <= counter_isco_multiple_rings_count <= 9:
                                        print(f"{Fore.LIGHTCYAN_EX} select_isco_orbital says:{Fore.LIGHTRED_EX} Incorrect Counter-ISCO Rings Count! must be between '1' and '9'!")
                                        print()
                                        continue
                                    else:
                                        print(f"{Fore.LIGHTGREEN_EX} {counter_isco_multiple_rings_count} Counter-ISCO Rings will be shown! {Fore.RESET}")
                                        print()
                                    break
                                except ValueError:
                                    print(f"{Fore.LIGHTCYAN_EX} select_isco_orbital says:{Fore.LIGHTRED_EX} Incorrect Counter-ISCO Rings Count! must be between '1' and '9'!")
                                    print()
                                    continue

                        case "3":
                            counter_isco_sphere = True
                            print(f"{Fore.LIGHTGREEN_EX} Counter-ISCO Sphere will be shown as Full Spherical! {Fore.RESET}")
                            print()
                    exit_first_loop = 1
                    break
                if exit_first_loop == 1:
                    break      
    return (
        co_isco_ring, 
        co_isco_ring_angle, 
        co_isco_multiple_rings, 
        co_isco_multiple_rings_count, 
        co_isco_sphere,
        counter_isco_ring, 
        counter_isco_ring_angle, 
        counter_isco_multiple_rings, 
        counter_isco_multiple_rings_count, 
        counter_isco_sphere
    )    

def isco_orbitals_info_print(bh_type, bh_mass, a_star, co_isco_ring, co_isco_ring_angle, co_isco_multiple_rings, co_isco_multiple_rings_count, co_isco_sphere, counter_isco_ring, counter_isco_ring_angle, counter_isco_multiple_rings, counter_isco_multiple_rings_count, counter_isco_sphere):
    print()
    print(Fore.BLUE + "==================================================" + Fore.RESET)
    print(Fore.RED + "Your Black Hole ISCO Orbitals" + Fore.RESET)

    match bh_type:
        case "schw" | "1":
            #single isco ring
            if co_isco_ring == True:
                print(Fore.LIGHTMAGENTA_EX + f"Single ISCO Ring:", end="\n")

                isco_sphere_radius(bh_type, bh_mass, a_star, "", np.radians(co_isco_ring_angle), True)

            #multiple isco rings
            elif co_isco_multiple_rings == True:
                print(Fore.LIGHTMAGENTA_EX + f"Multiple ISCO Rings:", end="\n")

                additional_angle = 180 / co_isco_multiple_rings_count
                angle = 0
                for i in range(co_isco_multiple_rings_count):
                    isco_sphere_radius(bh_type, bh_mass, a_star, "", np.radians(angle), True)
                    angle += additional_angle

            #isco sphere
            elif co_isco_sphere == True:
                print(Fore.LIGHTMAGENTA_EX + f"ISCO Sphere:", end="\n")

                additional_angle = 180 / 9
                angle = 0
                for i in range(9):
                    isco_sphere_radius(bh_type, bh_mass, a_star, "", np.radians(angle), True)
                    angle += additional_angle

        case "kerr" | "2":
            #co-rotating single isco ring
            if co_isco_ring == True:
                print(Fore.LIGHTMAGENTA_EX + f"Co-Rotating Single ISCO Ring:", end="\n")

                isco_sphere_radius(bh_type, bh_mass, a_star, "prograde", np.radians(co_isco_ring_angle), True)

            #co-rotating multiple isco rings
            elif co_isco_multiple_rings == True:
                print(Fore.LIGHTMAGENTA_EX + f"Co-Rotating Multiple ISCO Rings:", end="\n")

                additional_angle = 180 / co_isco_multiple_rings_count
                angle = 0
                for i in range(co_isco_multiple_rings_count):
                    isco_sphere_radius(bh_type, bh_mass, a_star, "prograde", np.radians(angle), True)
                    angle += additional_angle

            #co-rotating isco sphere
            elif co_isco_sphere == True:
                print(Fore.LIGHTMAGENTA_EX + f"Co-Rotating ISCO Sphere:", end="\n")

                additional_angle = 180 / 9
                angle = 0
                for i in range(9):
                    isco_sphere_radius(bh_type, bh_mass, a_star, "prograde", np.radians(angle), True)
                    angle += additional_angle

            #counter-rotating single isco ring
            if counter_isco_ring == True:
                print(Fore.LIGHTMAGENTA_EX + f"Counter-Rotating Single ISCO Ring:", end="\n")

                isco_sphere_radius(bh_type, bh_mass, a_star, "retrograde", np.radians(counter_isco_ring_angle), True)

            #counter-rotating isco rings
            if counter_isco_multiple_rings == True:

                print(Fore.LIGHTMAGENTA_EX + f"Counter-Rotating Multiple ISCO Rings:", end="\n")

                additional_angle = 180 / counter_isco_multiple_rings_count
                angle = 0
                for i in range(counter_isco_multiple_rings_count):
                    isco_sphere_radius(bh_type, bh_mass, a_star, "retrograde", np.radians(angle), True)
                    angle += additional_angle

            #co-rotating isco sphere
            elif counter_isco_sphere == True:
                print(Fore.LIGHTMAGENTA_EX + f"Counter-Rotating ISCO Sphere:", end="\n")

                additional_angle = 180 / 9
                angle = 0
                for i in range(9):
                    isco_sphere_radius(bh_type, bh_mass, a_star, "retrograde", np.radians(angle), True)
                    angle += additional_angle

def isco_orbitals_construct(bh_type, bh_mass, a_star, co_isco_ring, co_isco_ring_angle, co_isco_multiple_rings, co_isco_multiple_rings_count, co_isco_sphere, counter_isco_ring, counter_isco_ring_angle, counter_isco_multiple_rings, counter_isco_multiple_rings_count, counter_isco_sphere):
    match bh_type:
        case "schw" | "1":
            #single isco ring
            if co_isco_ring == True:
                isco_ring_draw(bh_type, bh_mass, a_star, "", np.radians(co_isco_ring_angle))
            
            #multiple isco rings
            elif co_isco_multiple_rings == True:
                additional_angle = 180 / co_isco_multiple_rings_count
                angle = 0
                for i in range(co_isco_multiple_rings_count):
                    isco_ring_draw(bh_type, bh_mass, a_star, "", np.radians(angle))
                    angle += additional_angle

            #isco sphere
            elif co_isco_sphere == True:
                isco_sphere_draw(bh_type, bh_mass, a_star, "")
                
        case "kerr" | "2":
            #co-rotating single isco ring
            if co_isco_ring == True:
                isco_ring_draw(bh_type, bh_mass, a_star, "prograde", np.radians(co_isco_ring_angle))

            #co-rotating multiple isco rings
            elif co_isco_multiple_rings == True:
                additional_angle = 180 / co_isco_multiple_rings_count
                angle = 0
                for i in range(co_isco_multiple_rings_count):
                    isco_ring_draw(bh_type, bh_mass, a_star, "prograde", np.radians(angle))
                    angle += additional_angle

            #co-rotating isco sphere
            elif co_isco_sphere == True:
                isco_sphere_draw(bh_type, bh_mass, a_star, "prograde")

            #counter-rotating single isco ring
            if counter_isco_ring == True:
                isco_ring_draw(bh_type, bh_mass, a_star, "retrograde", np.radians(counter_isco_ring_angle))


            #counter-rotating single isco rings
            elif counter_isco_multiple_rings == True:
                additional_angle = 180 / counter_isco_multiple_rings_count
                angle = 0
                for i in range(counter_isco_multiple_rings_count):
                    isco_ring_draw(bh_type, bh_mass, a_star, "retrograde", np.radians(angle))
                    angle += additional_angle

            #counter-rotating isco sphere
            elif counter_isco_sphere == True:
                isco_sphere_draw(bh_type, bh_mass, a_star, "retrograde")

#selecting object experience to show
def select_object_experience_type(bh_type, bh_mass, a_star):
    distance, direction, s_polar_angle, s_azimuthal_angle, velocity, v_polar_angle, v_azimuthal_angle = None, None, None, None, None, None, None
    while True:
        xp_type = input(Fore.LIGHTYELLOW_EX + "Choose the type of experience around this Monster!\n"
        + Fore.LIGHTBLACK_EX + "Enter '0' as no experience\n"
        "Enter '1' as Stationary\n"
        "Enter '2' as Orbiting\n"
        "Enter '3' as Free-Fall\n"
        "Enter '4' as Velocity-Initialized-Fall\n"
        + Fore.RESET + "==============================\n"
        ">>>")

        if xp_type not in ["0", "1", "2", "3", "4"]:
            print(f"{Fore.LIGHTCYAN_EX} select_experience_type says:{Fore.LIGHTRED_EX} Incorrect entry! must be in [0, 1, 2, 3, 4]")
            print()
            continue
        else:
            break

    match xp_type:
        case "0":
            pass
        case "1":
            s_polar_angle, s_azimuthal_angle = get_start_coordinates("1")
            print(Fore.RED + "----------------------------------------" + Fore.RESET)
            distance = get_distance("1", bh_type, bh_mass, a_star, s_polar_angle, direction)
            print(Fore.RED + "----------------------------------------" + Fore.RESET)
        case "2":
            s_polar_angle, s_azimuthal_angle = get_start_coordinates("2")
            print(Fore.RED + "----------------------------------------" + Fore.RESET)
            if bh_type in ["kerr", "2"]:
                direction = get_direction()
            print(Fore.RED + "----------------------------------------" + Fore.RESET)
            distance = get_distance("2", bh_type, bh_mass, a_star, s_polar_angle, direction)
            print(Fore.RED + "----------------------------------------" + Fore.RESET)
        case "3":
            s_polar_angle, s_azimuthal_angle = get_start_coordinates("3")
            print(Fore.RED + "----------------------------------------" + Fore.RESET)
            distance = get_distance("3", bh_type, bh_mass, a_star, s_polar_angle, direction)
            print(Fore.RED + "----------------------------------------" + Fore.RESET)
        case "4":
            s_polar_angle, s_azimuthal_angle = get_start_coordinates("4")
            print(Fore.RED + "----------------------------------------" + Fore.RESET)
            distance = get_distance("4", bh_type, bh_mass, a_star, s_polar_angle, direction)
            print(Fore.RED + "----------------------------------------" + Fore.RESET)
            velocity = get_velocity()
            print(Fore.RED + "----------------------------------------" + Fore.RESET)
            v_polar_angle, v_azimuthal_angle = get_velocity_coordinates()
            print(Fore.RED + "----------------------------------------" + Fore.RESET)

    return (
        xp_type,
        distance,
        direction, 
        s_polar_angle, 
        s_azimuthal_angle, 
        velocity, 
        v_polar_angle, 
        v_azimuthal_angle
    )

def get_start_coordinates(xp_type):
    match xp_type:
        case "1":
            first_polar_message = "in which Polar Angle you want to Stay near this Monster?"
            first_azimuthal_message = "in which Azimuthal Angle you want to Stay near this Monster?"
            second_polar_message = "Stationary Object at the Polar Angle of:"
            second_azimuthal_message = "Stationary Object at the Azimuthal Angle of:"
        case "2":
            first_polar_message = "in which Polar Angle you want to Orbit near this Monster?"
            first_azimuthal_message = "in which Azimuthal Angle you want to Orbit near this Monster?"
            second_polar_message = "Orbiting Object at the Polar Angle of:"
            second_azimuthal_message = "Orbiting Object at the Azimuthal Angle of:"
        case "3" | "4":
            first_polar_message = "in which Polar Angle you want to Fall Toward this Monster?"
            first_azimuthal_message = "in which Azimuthal Angle you want to Fall Toward this Monster?"
            second_polar_message = "Falling Object at the Polar Angle of:"
            second_azimuthal_message = "Falling Object at the Azimuthal Angle of:"
    #get polar angle  
    while True:
        try:
            polar_angle = int(input(Fore.LIGHTYELLOW_EX + f"{first_polar_message}\n"
            + Fore.LIGHTBLACK_EX + "Must be between 0° and 180°\n"
            + Fore.RESET + "==============================\n>>>"))
            if not 0 <= polar_angle <= 180:
                print(f"{Fore.LIGHTCYAN_EX} get_start_coordinates says:{Fore.LIGHTRED_EX} Incorrect Polar Angle! must be between '0' and '180'!")
                print()
                continue

            break

        except ValueError:
            print(f"{Fore.LIGHTCYAN_EX} get_start_coordinates says:{Fore.LIGHTRED_EX} Incorrect Entry! must be number!")
            print()
            continue

    #get azimuthal angle  
    while True:
        try:
            azimuthal_angle = int(input(Fore.LIGHTYELLOW_EX + f"{first_azimuthal_message}\n"
            + Fore.LIGHTBLACK_EX + "Must be between 0° and 360°\n" 
            + Fore.RESET + "==============================\n>>>"))
            if not 0 <= azimuthal_angle <= 360:
                print(f"{Fore.LIGHTCYAN_EX} get_start_coordinates says:{Fore.LIGHTRED_EX} Incorrect Azimuthal Angle! must be between '0' and '360'!")
                print()
                continue

            break

        except ValueError:
            print(f"{Fore.LIGHTCYAN_EX} get_start_coordinates says:{Fore.LIGHTRED_EX} Incorrect Entry! must be number!")
            print()
            continue

    print(f"{Fore.LIGHTGREEN_EX} {second_polar_message} {Fore.RESET}{polar_angle}°.")
    print(f"{Fore.LIGHTGREEN_EX} {second_azimuthal_message} {Fore.RESET}{azimuthal_angle}°.")

    return np.radians(polar_angle), np.radians(azimuthal_angle)

def get_direction():
    while True:
        direction = input(Fore.LIGHTYELLOW_EX + f"in which Direction you want to orbit this Monster?\n"
        + Fore.LIGHTBLACK_EX + "Enter '1' as Prograde\n"
        "Enter '2' as Retrograde\n"
        + Fore.RESET + "==============================\n"
        ">>>")

        if direction not in ["1", "2"]:
            print(f"{Fore.LIGHTCYAN_EX} get_experience_type says:{Fore.LIGHTRED_EX} Incorrect Direction! must be '1' or '2'!")
            print()
            continue
        else:
            print(f"{Fore.LIGHTGREEN_EX} Prograde Orbting Object. {Fore.RESET}") if direction == "1" else print(f"{Fore.LIGHTGREEN_EX} Retrograde Orbiting Object. {Fore.RESET}")
            print()
        break
    
    return "prograde" if direction == "1" else "retrograde"

def get_distance(xp_type, bh_type, bh_mass, a_star, s_polar_angle, direction):
    match xp_type:
        case "1":
            first_message = "How near you want to be at the Event Horizon?"
            second_message = "Stationary Object at the Boyer-Lindquist Distance of "
            #calculating minimum safe place outside the ergosphere
            safe_distance = outer_ergosphere_radius(bh_mass, a_star, s_polar_angle, False) - outer_event_horizon_radius(bh_type, bh_mass, a_star, s_polar_angle, False)
            hint_message = (
                f"{Fore.LIGHTYELLOW_EX}Hint : Minimum Safe Boyer-Lindquist Distance from Event Horizon at a Polar Angle of {Fore.LIGHTRED_EX}{int(np.round(np.degrees(s_polar_angle)))}°{Fore.LIGHTYELLOW_EX} should be more than: {Fore.RESET}{length_unit_manager(safe_distance)}.\n"
                f"{Fore.LIGHTYELLOW_EX} Cartesian Equivalent is {Fore.RESET}{length_unit_manager(cartesian_radius(bh_mass, outer_ergosphere_radius(bh_mass, a_star, s_polar_angle, False), a_star, s_polar_angle) - cartesian_radius(bh_mass, outer_event_horizon_radius(bh_type, bh_mass, a_star, s_polar_angle, False), a_star, s_polar_angle))}\n"
                f"{Fore.LIGHTYELLOW_EX}Distance should be greater than the Event Horizon Radius.{Fore.RESET}"
            ) 
        case "2":
            first_message = "How near you want to orbit around the Event Horizon?"
            second_message = "Orbiting Object at the Boyer-Lindquist Distance of"
            #calculating minimum safe place outside the photon ring
            safe_distance = photon_sphere_radius(bh_type, bh_mass, a_star, direction, s_polar_angle, False) - outer_event_horizon_radius(bh_type, bh_mass, a_star, s_polar_angle, False)
            hint_message = (
                f"{Fore.LIGHTYELLOW_EX}Hint : Minimum Safe Boyer-Lindquist Distance from Event Horizon at a Polar Angle of {Fore.LIGHTRED_EX}{int(np.round(np.degrees(s_polar_angle)))}°{Fore.LIGHTYELLOW_EX} should be more than: {Fore.RESET}{length_unit_manager(safe_distance)}.\n"
                f"{Fore.LIGHTYELLOW_EX} Cartesian Equivalent is {Fore.RESET}{length_unit_manager(cartesian_radius(bh_mass, photon_sphere_radius(bh_type, bh_mass, a_star, direction, s_polar_angle, False), a_star, s_polar_angle) - cartesian_radius(bh_mass, outer_event_horizon_radius(bh_type, bh_mass, a_star, s_polar_angle, False), a_star, s_polar_angle))}\n"
                f"{Fore.LIGHTYELLOW_EX}Distance should be greater than the Photon Ring Radius.{Fore.RESET}"
            )
        case "3":
            first_message = "From which distance you want to fall toward the Event Horizon?"
            second_message = "Falling Object at the Boyer-Lindquist Distance of"
        case "4":
            first_message = "From which distance you want to ignite toward the Event Horizon?"
            second_message = "Speeding Object at the Boyer-Lindquist Distance of"

    while True:
        distance = input(Fore.LIGHTYELLOW_EX + f"{first_message}\n"
        + Fore.LIGHTBLACK_EX + "use 'm' as Meter\n"
        "use 'km' as KiloMeter\n"
        "use 'au' as Astronomical Unit\n"
        "use 'ly' as Light Years\n"
        "or enter without any to use 'm' as default\n"
        f"{hint_message}\n"
        + Fore.RESET + "==============================\n"
        ">>>")
        try:
            distance_type = ""
            for i in distance:
                if i.isalpha():
                    distance_type = distance_type + str(i)
                    distance = distance.replace(i, "")

            #check unit
            if distance_type == "": distance_type = "m"
            if distance_type not in distance_type_list:
                print(f" {Fore.LIGHTCYAN_EX} get-distance says: {Fore.LIGHTRED_EX} unknown unit type: {Fore.RESET} {distance_type}")
                continue

            if float(distance) < 0:
                print(f" {Fore.LIGHTCYAN_EX} get-distance says: {Fore.LIGHTRED_EX} you CAN NOT be inside the Event Horizon!! {Fore.RESET}")
                distance_type = ""
                print()
                continue

            elif float(distance) == 0:
                print(f" {Fore.LIGHTCYAN_EX} get-distance says: {Fore.LIGHTYELLOW_EX} Warning! you're at the edge of the Event Horizon!!! {Fore.RESET}")                  

            print(f"{Fore.LIGHTGREEN_EX} {second_message}{Fore.RESET} {distance}{distance_type} {Fore.LIGHTGREEN_EX} from the Event Horizon! {Fore.RESET}")
            distance = (float(distance) * distance_type_list[distance_type]) + outer_event_horizon_radius(bh_type, bh_mass, a_star, "", False)
            print(f"{Fore.LIGHTGREEN_EX} Boyer-Lindquist Distance from the Black Hole's Center: {Fore.RESET}{length_unit_manager(distance)}.")
            break

        except ValueError:
            print(Fore.LIGHTRED_EX + "distance must have number!" + Fore.RESET)
            distance_type = ""
            print()
            continue

    return float(distance)

def get_velocity():
    while True:
        velocity = input(Fore.LIGHTYELLOW_EX + f"Enter your initial velocity.\n"
        + Fore.LIGHTBLACK_EX + "use 'cm/s' as CentiMeter/s\n"
        "use 'm/s' as Meter/s\n"
        "use 'sos/s' as Speed of Sound/s\n"
        "use 'km/s' as KiloMeter/s\n"
        "use 'c/s' as Speed of Light/s\n"
        "or enter without any to use 'm/s' as default\n"
        + Fore.RESET + "==============================\n"
        ">>>")
        try:
            velocity_type = ""
            for i in velocity:
                if i.isalpha():
                    velocity_type = velocity_type + str(i)
                    velocity = velocity.replace(i, "")

            #check unit
            if velocity_type == "": velocity_type = "m/s"
            if velocity_type not in velocity_type_list:
                print(f" {Fore.LIGHTCYAN_EX} get_velocity says: {Fore.LIGHTRED_EX} unknown unit type: {Fore.RESET} {velocity_type}")
                continue

            if float(velocity) <= 0:
                print(f" {Fore.LIGHTCYAN_EX} get_velocity says: {Fore.LIGHTRED_EX} speed can not be Zero or Negative!! {Fore.RESET}")
                velocity_type = ""
                print()
                continue               

            print(f"{Fore.LIGHTGREEN_EX} Your Initial Velocity: {velocity}{velocity_type}! {Fore.RESET}")
            print()
            break

        except ValueError:
            print(Fore.LIGHTRED_EX + "Velocity must have number!" + Fore.RESET)
            velocity_type = ""
            print()
            continue

    return float(velocity) * velocity_type_list[velocity_type]

def get_velocity_coordinates():
    #get polar angle  
    while True:
        try:
            polar_angle = int(input(Fore.LIGHTYELLOW_EX + f"in which Polar Angle you want to Initiate your engines?\n"
            + Fore.LIGHTBLACK_EX + "Must be between 0° and 180°\n"))
            if not 0 <= polar_angle <= 180:
                print(f"{Fore.LIGHTCYAN_EX} get_velocity_coordinates says:{Fore.LIGHTRED_EX} Incorrect Polar Angle! must be between '0' and '180'!")
                print()
                continue
            else:
                print(f"{Fore.LIGHTGREEN_EX} you want to start your engines at the Polar Angle of {polar_angle}°! {Fore.RESET}")
                print()
            break

        except ValueError:
            print(f"{Fore.LIGHTCYAN_EX} get_velocity_coordinates says:{Fore.LIGHTRED_EX} Incorrect Entry! must be number!")
            print()
            continue

    #get azimuthal angle  
    while True:
        try:
            azimuthal_angle = int(input(Fore.LIGHTYELLOW_EX + f"in which Azimuthal Angle you want to Initiate your engines?\n"
            + Fore.LIGHTBLACK_EX + "Must be between 0° and 360°\n"))
            if not 0 <= azimuthal_angle <= 360:
                print(f"{Fore.LIGHTCYAN_EX} get_velocity_coordinates says:{Fore.LIGHTRED_EX} Incorrect Azimuthal Angle! must be between '0' and '360'!")
                print()
                continue
            else:
                print(f"{Fore.LIGHTGREEN_EX} you want to start your engines at the Azimuthal Angle of {azimuthal_angle}°! {Fore.RESET}")
                print()
            break

        except ValueError:
            print(f"{Fore.LIGHTCYAN_EX} get_velocity_coordinates says:{Fore.LIGHTRED_EX} Incorrect Entry! must be number!")
            print()
            continue

    return np.radians(polar_angle), np.radians(azimuthal_angle)

def object_experience(bh_type, bh_mass, a_star, xp_type, distance, direction, s_polar_angle, s_azimuthal_angle, velocity, v_polar_angle, v_azimuthal_angle):
    print()
    print(Fore.BLUE + "==================================================" + Fore.RESET)

    match xp_type:
        case "1":
            print(Fore.RED + "Stationary Object Statistics" + Fore.RESET)
            stationary_object_xp_info_print(bh_type, bh_mass, a_star, distance, s_polar_angle, s_azimuthal_angle)
            stationary_object_xp_draw(bh_type, bh_mass, a_star, distance, s_polar_angle, s_azimuthal_angle)
        case "2":
            print(Fore.RED + "Orbiting Object Statistics" + Fore.RESET)
            coordinate_data = stationary_object_xp_draw(bh_type, bh_mass, a_star, distance, s_polar_angle, s_azimuthal_angle)
            orbiting_object_xp_draw(bh_type, bh_mass, a_star, distance, direction, s_polar_angle, coordinate_data[0], coordinate_data[1], coordinate_data[2], coordinate_data[3], coordinate_data[4], coordinate_data[5])
        case "3":
            ...
        case "4":
            ...

#converters
def length_unit_manager(item):
    for limit, function, unit in length_type_list:
        if item < limit:
            length = function(item)
            return f"{length:,.5f} {unit}"
              
def speed_unit_manager(velocity):
    for limit, function, unit in speed_type_list:
        if velocity < limit:
            value = function(velocity)
            return f"{value:,.9f} {unit}"
    
def cartesian(bh_mass, radius, theta, phi, a_star):
    a = a_star * gravitational_radius(bh_mass)

    x = np.sqrt((radius ** 2) + (a ** 2)) * np.sin(theta) * np.cos(phi)
    y = np.sqrt((radius ** 2) + (a ** 2)) * np.sin(theta) * np.sin(phi)
    z = radius * np.cos(theta)

    return x, y, z

def cartesian_radius(bh_mass, radius, a_star, theta):
    r_cartesian = np.sqrt((radius ** 2) + (((a_star) * (gravitational_radius(bh_mass))) ** 2) * (np.sin(theta) ** 2))
    return r_cartesian

#equations
def photon_equation(bh_mass, a_star, radius, direction, i):
    match direction:
        case "prograde":
            i = (np.pi / 2) - i
        case "retrograde":
            i = (np.pi / 2) + i

    a = a_star * gravitational_radius(bh_mass)
    m = gravitational_radius(bh_mass)

    if 0 < a_star <= 1:
        xi = (((radius ** 2) * (radius - (3 * m))) + ((a ** 2) * (radius + m))) / (a * (m - radius))
        eta = ((radius ** 3) * ((4 * (a ** 2) * m) - radius * (radius - (3 * m)) ** 2)) / ((a ** 2) * (m - radius) ** 2)
            
        equation = np.cos(i) - (xi) / np.sqrt((xi ** 2) + (eta))

    return equation

def isco_equation(bh_mass, a_star, radius, direction, i):
    match direction:
        case "prograde":
            i = (np.pi / 2) - i
        case "retrograde":
            i = (np.pi / 2) + i

    a = a_star * gravitational_radius(bh_mass)
    m = gravitational_radius(bh_mass)

    #Z(r)
    zr = ((radius * (radius - 6 * m)) ** 2) - ((a ** 2) * (((2 * radius) * ((3 * radius) + (14 * m))) - (9 * (a ** 2))))
    #X(r)
    xr = ((a ** 2) * (((a ** 2) * ((3 * (a ** 2)) + ((4 * radius) * ((2 * radius) - (3 * m))))) + ((radius ** 2) * (((15 * radius) * (radius - (4 * m))) + (28 * (m ** 2)))))) - ((6 * (radius ** 4)) * ((radius ** 2) - (4 * (m ** 2))))
    #Y(r)
    yr = ((a ** 4) * ((a ** 4) + ((radius ** 2) * (((7 * radius) * ((3 * radius) - (4 * m))) + (36 * (m ** 2)))))) + (((6 * radius) * (radius - (2 * m))) * ((a ** 6) + ((2 * (radius ** 3)) * (((a ** 2) * ((3 * radius) + (2 * m))) + ((3 * (radius ** 2)) * (radius - (2 * m)))))))

    #S(r)
    equation = ((radius ** 8) * zr) + (((a ** 2) * (1 - (np.cos(i) ** 2))) * (((a ** 2) * (1 - (np.cos(i) ** 2)) * yr) - ((2 * (radius ** 4)) * xr)))

    return equation

def polar_isco_equation(bh_mass, a_star, radius):
    a = a_star * gravitational_radius(bh_mass)
    m = gravitational_radius(bh_mass)

    #P(r)
    equation = ((radius ** 3) * (((radius ** 2) * (radius - (6 * m))) + ((a ** 2) * ((3 * radius) + (4 * m))))) + ((a ** 4) * (((3 * radius) * (radius - (2 * m))) + (a ** 2)))

    return equation

#black hole main structure
def gravitational_radius(bh_mass):
    gravitational_radius = (bh_mass * gravitational_constant) / c2
    return gravitational_radius

def singularity_radius(bh_mass, a_star, show):
    match a_star:
        case 0:
            singularity_radius = 0
            if show == True:
                print(Fore.LIGHTYELLOW_EX + f"Point Singularity at the Very Center with no Radius! {Fore.RESET}", end="\n")
                print(Fore.RED + "----------------------------------------" + Fore.RESET)  
        case _:
            singularity_radius = a_star * gravitational_radius(bh_mass)
            if show == True:
                print(Fore.LIGHTYELLOW_EX + f"Ring Singularity Boyer-Lindquist Radius: {Fore.RESET} 0!", end="\n")
                print("         ********************            ")
                print(Fore.LIGHTYELLOW_EX + f"Ring Singularity Cartesian Radius: {Fore.RESET} {length_unit_manager(singularity_radius)} \n {Fore.LIGHTYELLOW_EX} Diameter: {Fore.RESET} {length_unit_manager(singularity_radius * 2)}!" + Fore.RESET)  
                print(Fore.RED + "----------------------------------------" + Fore.RESET) 
    return singularity_radius

def singularity_draw(bh_mass, a_star, singularity_radius, outer_event_horizon_radius):
    match a_star:
        case 0:
            #from north pole to south pole of the singularity
            theta = np.linspace(0, np.pi, 100)

            #a turn around the point singularity
            phi = np.linspace(0, 2 * np.pi, 100)

            #create dots
            theta, phi = np.meshgrid(theta, phi) 

            #calling cartesian function to convert point singularity coordinates from boyer lindquist to cartesian
            x, y, z = cartesian(bh_mass, outer_event_horizon_radius / 70, theta, phi, a_star)

            text = "The Point Singularity Radius:"

            #creatig the black hole's point singularity
            ax.plot_surface(x, y, z, color="yellow", alpha=0.7)

            #drawing the arrow for poit singularity
            ax.plot([0, 0], [0, -outer_event_horizon_radius * 1.3], [0, 0], color="yellow")

            #adding the text for poit singularity
            ax.text(0, -outer_event_horizon_radius * 1.3, 0, f"{text} " f"{length_unit_manager(singularity_radius)}", color="yellow", path_effects=[pe.withStroke(linewidth=2, foreground="black")], fontsize=10)

            #setting base for the coordinates of the equatorial informations
            limit = outer_event_horizon_radius / 70

            ax.set_xlim(-limit * 100, limit * 100)
            ax.set_ylim(-limit * 100, limit * 100)
            ax.set_zlim(-limit * 100, limit * 100)
        case _:
            # The ring singularity lies in the equatorial plane
            theta = np.pi / 2       

            #a turn around the ring singularity
            phi = np.linspace(0, 2 * np.pi, 100)

            #create the mesh of dots for the surface
            theta, phi = np.meshgrid(theta, phi) 

            #calling ring_singularity_cartesian function to convert ring singularity coordinates from boyer lindquist to cartesian
            x, y, z = cartesian(bh_mass, 0, theta, phi, a_star)

            text = "The Ring Singularity Cartesian Radius:"

            #creatig the black hole's ring singularity
            ax.plot(x, y, z, color="yellow", linewidth=5)

            angle = np.radians(210)

            x = singularity_radius * np.sin(angle)
            y = singularity_radius * np.cos(angle)
            z = 0   

            #drawing the arrow for ring singularity
            ax.plot([x, x + outer_event_horizon_radius * 1.8 * np.sin(angle)], [y, y + outer_event_horizon_radius * 1.8 * np.cos(angle)], [z, z], color="yellow")

            #adding the text for ring singularity
            ax.text(x + outer_event_horizon_radius * 1.8 * np.sin(angle), y + outer_event_horizon_radius * 1.8 * np.cos(angle), z, f"{text} " f"{length_unit_manager(singularity_radius)}", color="yellow", path_effects=[pe.withStroke(linewidth=2, foreground="black")], fontsize=10)

            limit = singularity_radius * 1.5

            ax.set_xlim(-limit * 4, limit * 4)
            ax.set_ylim(-limit * 4, limit * 4)
            ax.set_zlim(-limit * 4, limit * 4)

    ax.set_box_aspect([1, 1, 1])

    #show the drawings
    plt.pause(0.5)

def inner_ergosphere_radius(bh_mass, a_star, theta, show):
    inner_ergosphere_radius = gravitational_radius(bh_mass) * (1 - np.sqrt(1 - ((a_star ** 2) * (np.cos(theta) ** 2))))
    if show == True:
        print(Fore.LIGHTGREEN_EX + f"Inner Ergosphere Boyer-Linquist Radius: {Fore.RESET} {length_unit_manager(inner_ergosphere_radius)}", end="\n")
        print(Fore.LIGHTGREEN_EX + f"Inner Ergosphere Boyer-Linquist Diameter: {Fore.RESET} {length_unit_manager(inner_ergosphere_radius * 2)}", end="\n")
        print("         ********************            ")
        print(Fore.LIGHTCYAN_EX + f"Inner Ergosphere Cartesian Radius: {Fore.RESET} {length_unit_manager(cartesian_radius(bh_mass, inner_ergosphere_radius, a_star, theta))}", end="\n")
        print(Fore.LIGHTCYAN_EX + f"Inner Ergosphere Cartesian Diameter: {Fore.RESET} {length_unit_manager(cartesian_radius(bh_mass, inner_ergosphere_radius, a_star, theta) * 2)}", end="\n")
        print(Fore.RED + "----------------------------------------" + Fore.RESET) 
    return inner_ergosphere_radius

def inner_ergosphere_draw(bh_mass, a_star):
    if a_star != 0:
        #from north pole to south pole of the inner ergosphere
        theta = np.linspace(0, np.pi, 100)

        #a turn around the inner ergosphere
        phi = np.linspace(0, 2 * np.pi, 100)

        #create dots
        theta, phi = np.meshgrid(theta, phi)

        #creating an array of radius based on each theta
        r = inner_ergosphere_radius(bh_mass, a_star, theta, False)

        x, y, z = cartesian(bh_mass, r, theta, phi, a_star)

        ax.plot_surface(x, y, z, color="#460000", alpha=0.2, linewidth=3)

        #setting base for the coordinates of the polar informations
        polar_limit = inner_ergosphere_radius(bh_mass, a_star, np.radians(0), False)

        #setting base for the coordinates of the equatorial informations
        equatorial_limit = cartesian_radius(bh_mass, inner_ergosphere_radius(bh_mass, a_star, np.radians(90), False), a_star, np.radians(90))

        angle = np.radians(270)

        #setting coordinates for text and arrow for inner ergosphere cartesian polar radius
        x = polar_limit * np.sin(angle)
        y = polar_limit * np.cos(angle)
        z = polar_limit

        #drawing an arrow for inner ergosphere cartesian polar radius
        ax.plot([0, polar_limit], [0, 0], [z, z * 3], color="#460000")

        #adding an text for inner ergosphere cartesian polar radius
        ax.text(polar_limit, 0, z * 3, f"Inner Ergosphere Cartesian Polar Radius: " f"{length_unit_manager(polar_limit)}", path_effects=[pe.withStroke(linewidth=2, foreground="white")], color="#460000",fontweight="bold", fontsize=10)


        #setting coordinates for text and arrow for inner ergosphere cartesian equatorial radius
        x = equatorial_limit * np.sin(angle)
        y = equatorial_limit * np.cos(angle)
        z = 0

        #drawing an arrow for inner ergosphere cartesian equatorial radius
        ax.plot([x, x * 3], [y, y * 3], [z, z], color="#460000")

        #adding the text for inner ergosphere cartesian equatorial radius
        ax.text(x * 3, y * 3, z, f"Inner Ergosphere Cartesian Equatorial Radius: " f"{length_unit_manager(equatorial_limit)}", path_effects=[pe.withStroke(linewidth=2, foreground="white")], color="#460000",fontweight="bold", fontsize=10)


        ax.set_xlim(-equatorial_limit * 1.2, equatorial_limit * 1.2)
        ax.set_ylim(-equatorial_limit * 1.2, equatorial_limit * 1.2)
        ax.set_zlim(-equatorial_limit * 1.2, equatorial_limit * 1.2)

        ax.set_box_aspect([1, 1, 1])

        #show the drawings
        plt.pause(0.5)

def inner_event_horizon_radius(bh_type, bh_mass, a_star, theta, show):
    if bh_type == "kerr" or bh_type == "2":
        inner_radius = gravitational_radius(bh_mass) * (1 - sqrt(1 - pow(a_star,2)))
        if show == True:
            print(Fore.LIGHTGREEN_EX + f"Inner Event Horizon Boyer-Lindquist Radius: {Fore.RESET} {length_unit_manager(inner_radius)}", end="\n")
            print(Fore.LIGHTGREEN_EX + f"Inner Event Horizon Boyer-Lindquist Diameter: {Fore.RESET} {length_unit_manager(inner_radius * 2)}", end="\n")
            print("         ********************            ")
            print(Fore.LIGHTCYAN_EX + f"Inner Event Horizon Cartesian Radius: {Fore.RESET} {length_unit_manager(cartesian_radius(bh_mass, inner_radius, a_star, theta))}", end="\n")
            print(Fore.LIGHTCYAN_EX + f"Inner Event Horizon Cartesian Diameter: {Fore.RESET} {length_unit_manager(cartesian_radius(bh_mass, inner_radius, a_star, theta) * 2)}", end="\n")
            print(Fore.RED + "----------------------------------------" + Fore.RESET)  

        return inner_radius

def inner_event_horizon_draw(bh_type, bh_mass, a_star, inner_horizon_radius):
    #from north pole to south pole of the inner horizon
    theta = np.linspace(0, np.pi, 100)

    #a turn around the inner horizon
    phi = np.linspace(0, 2 * np.pi, 100)

    #create dots
    theta, phi = np.meshgrid(theta, phi)   

    #calling cartesian function to convert inner horizon coordinates from boyer-lindquist to cartesian
    x, y, z = cartesian(bh_mass, inner_horizon_radius, theta, phi, a_star)

    #creatig the black hole's inner event horizon
    ax.plot_surface(x, y, z, color="black", alpha=0.4)

    #setting base for the coordinates of the polar informations
    polar_limit = inner_event_horizon_radius(bh_type, bh_mass, a_star, np.radians(0), False)

    #setting base for the coordinates of the equatorial informations
    equatorial_limit = cartesian_radius(bh_mass, inner_horizon_radius, a_star, np.radians(90))

    angle = np.radians(55)

    #setting coordinates for text and arrow for inner horizon cartesian polar radius
    x = polar_limit * np.tan(angle)
    z = polar_limit * np.cos(angle)

    #drawing an arrow for inner ergosphere cartesian polar radius
    ax.plot([0, x], [0, 0], [polar_limit, z * 4], color="black")

    #adding an text for inner ergosphere cartesian polar radius
    ax.text(x, 0, z * 4, f"Inner Event Horizon Cartesian Polar Radius: " f"{length_unit_manager(polar_limit)}", path_effects=[pe.withStroke(linewidth=2, foreground="white")], color="black", fontweight="bold", fontsize=10)

    #setting coordinates for text and arrow for inner horizon cartesian equatorial radius
    x = equatorial_limit * np.sin(angle)
    y = equatorial_limit * np.cos(angle)
    z = 0

    #drawing an arrow for inner horizon cartesian equatorial radius
    ax.plot([x, x * 1.7], [y, y * 1.7], [z ,z], color="black")

    #adding the text for inner horizon cartesian equatorial radius
    ax.text(x * 1.7, y * 1.7, z, f"Inner Event Horizon Cartesian Equatorial Radius: " f"{length_unit_manager(cartesian_radius(bh_mass, inner_horizon_radius, a_star, np.radians(90)))}", path_effects=[pe.withStroke(linewidth=2, foreground="white")], color="black", fontweight="bold", fontsize=10)

    ax.set_xlim(-1 * equatorial_limit * 3, equatorial_limit * 3)
    ax.set_ylim(-1 * equatorial_limit * 3, equatorial_limit * 3)
    ax.set_zlim(-1 * equatorial_limit * 3, equatorial_limit * 3)

    ax.set_box_aspect([1, 1, 1])

    #show the drawings
    plt.pause(0.5)
    
def outer_event_horizon_radius(bh_type, bh_mass, a_star, theta, show):
    match bh_type:
        case "schw" | "1":
            outer_radius = ((2 * bh_mass * gravitational_constant) / c2)
            if show == True:
                print(Fore.LIGHTGREEN_EX + f"Event Horizon Radius: {Fore.RESET} {length_unit_manager(outer_radius)}", end="\n")
                print(Fore.LIGHTGREEN_EX + f"Event Horizon Diameter: {Fore.RESET} {length_unit_manager(outer_radius * 2)}", end="\n")  
                print(Fore.RED + "----------------------------------------" + Fore.RESET)

        case "kerr" | "2":  
            outer_radius = gravitational_radius(bh_mass) * (1 + sqrt(1 - pow(a_star,2)))
            if show == True:
                print(Fore.LIGHTGREEN_EX + f"Outer Event Horizon Boyer-Lindquist Radius: {Fore.RESET} {length_unit_manager(outer_radius)}", end="\n")
                print(Fore.LIGHTGREEN_EX + f"Outer Event Horizon Boyer-Lindquist Diameter: {Fore.RESET} {length_unit_manager(outer_radius * 2)}", end="\n")  
                print("         ********************            ")
                print(Fore.LIGHTCYAN_EX + f"Outer Event Horizon Cartesian Radius: {Fore.RESET} {length_unit_manager(cartesian_radius(bh_mass, outer_radius, a_star, theta))}", end="\n")
                print(Fore.LIGHTCYAN_EX + f"Outer Event Horizon Cartesian Diameter: {Fore.RESET} {length_unit_manager(cartesian_radius(bh_mass, outer_radius, a_star, theta) * 2)}", end="\n")
                print(Fore.RED + "----------------------------------------" + Fore.RESET)

    return outer_radius
    
def outer_event_horizon_draw(bh_mass, bh_type, a_star, outer_horizon_radius):
    #from north pole to south pole of the outer horizon
    theta = np.linspace(0, np.pi, 100)

    #a turn around the outer horizon
    phi = np.linspace(0, 2 * np.pi, 100)

    #create dots
    theta, phi = np.meshgrid(theta, phi)   

    #calling cartesian function to convert outer horizon coordinates from boyer-lindquist to cartesian
    x, y, z = cartesian(bh_mass, outer_horizon_radius, theta, phi, a_star)

    angle = np.radians(35)

    if bh_type == "schw" or bh_type == "1":
        #creatig the black hole's outer event horizon
        ax.plot_surface(x, y, z, color="black", alpha=0.7)

        #setting base for the coordinates of the equatorial informations
        equatorial_limit = outer_horizon_radius * 1.2

        #setting coordinates for text and arrow for outer horizon radius
        x = outer_horizon_radius * np.sin(angle)
        y = outer_horizon_radius * np.cos(angle)
        z = 0

        #drawing an arrow for the horizon radius
        ax.plot([x, x * 1.5], [y, y * 1.5], [z, z], color="darkgray")

        #adding the text for the horizon radius
        ax.text(x * 1.5, y * 1.5, z, f"Event Horizon Radius: " f"{length_unit_manager(outer_horizon_radius)}", path_effects=[pe.withStroke(linewidth=2, foreground="white")], color="black", fontweight="bold", fontsize=10)

    elif bh_type == "kerr" or bh_type == "2":
        #creatig the black hole's outer event horizon
        ax.plot_surface(x, y, z, color="black", alpha=0.5)

        #setting base for the coordinates of the polar informations
        polar_limit = outer_horizon_radius

        #setting base for the coordinates of the equatorial informations
        equatorial_limit = cartesian_radius(bh_mass, outer_horizon_radius, a_star, np.radians(90))

        #setting coordinates for text and arrow for outer horizon cartesian polar radius
        x = polar_limit * np.tan(angle)
        z = polar_limit * np.cos(angle)

        #drawing an arrow for inner ergosphere cartesian polar radius
        ax.plot([0, -x], [0, 0], [polar_limit, z * 1.7], color="black")

        #adding an text for inner ergosphere cartesian polar radius
        ax.text(-x, 0, z * 1.7, f"Outer Event Horizon Cartesian Polar Radius: " f"{length_unit_manager(polar_limit)}", path_effects=[pe.withStroke(linewidth=2, foreground="white")], color="black", fontweight="bold", fontsize=10)

        #setting coordinates for text and arrow for outer horizon cartesian equatorial radius
        x = equatorial_limit * np.sin(angle)
        y = equatorial_limit * np.cos(angle)
        z = 0   

        #drawing an arrow for outer horizon cartesian equatorial radius
        ax.plot([x, x * 1.4], [y, y * 1.4], [z, z], color="black")

        #adding the text for outer horizon cartesian equatorial radius
        ax.text(x * 1.4, y * 1.4, z, f"Outer Event Horizon Cartesian Equatorial Radius: " f"{length_unit_manager(cartesian_radius(bh_mass, outer_horizon_radius, a_star, np.radians(90)))}", path_effects=[pe.withStroke(linewidth=2, foreground="white")], color="black", fontweight="bold", fontsize=10)

    ax.set_xlim(-1 * equatorial_limit * 1.2, equatorial_limit * 1.2)
    ax.set_ylim(-1 * equatorial_limit * 1.2, equatorial_limit * 1.2)
    ax.set_zlim(-1 * equatorial_limit * 1.2, equatorial_limit * 1.2)

    ax.set_box_aspect([1, 1, 1])

    #show the drawings
    plt.pause(0.5)

def outer_ergosphere_radius(bh_mass, a_star, theta, show):
    outer_ergosphere_radius = gravitational_radius(bh_mass) * (1 + np.sqrt(1 - a_star ** 2 * np.cos(theta) ** 2))
    outer_ergosphere_radius_from_horizon = cartesian_radius(bh_mass, outer_ergosphere_radius, a_star, theta) - cartesian_radius(bh_mass, outer_event_horizon_radius("kerr", bh_mass, a_star, theta, False), a_star, theta)
    if show == True:
        print(Fore.LIGHTGREEN_EX + f"Outer Ergosphere Boyer-Lindquist Radius: {Fore.RESET} {length_unit_manager(outer_ergosphere_radius)}", end="\n")
        print(Fore.LIGHTGREEN_EX + f"Outer Ergosphere Boyer-Lindquist Diamater: {Fore.RESET} {length_unit_manager(outer_ergosphere_radius * 2)}", end="\n")
        print("         ********************            ")
        print(Fore.LIGHTCYAN_EX + f"Outer Ergosphere Cartesian Radius: {Fore.RESET} {length_unit_manager(cartesian_radius(bh_mass, outer_ergosphere_radius, a_star, theta))}", end="\n")
        print(Fore.LIGHTCYAN_EX + f"Outer Ergosphere Cartesian Diameter: {Fore.RESET} {length_unit_manager(cartesian_radius(bh_mass, outer_ergosphere_radius, a_star, theta) * 2)}", end="\n")
        print("         ********************            ")
        print(Fore.LIGHTRED_EX + f"Extended about: {Fore.RESET} {length_unit_manager(outer_ergosphere_radius_from_horizon)} {Fore.LIGHTRED_EX} from the Outer Event Horizon at the Equator.", end="\n")
        print(Fore.RED + "----------------------------------------" + Fore.RESET)

    return outer_ergosphere_radius

def outer_ergosphere_draw(bh_type, bh_mass, a_star):
    if bh_type == "kerr" or bh_type == "2":
        #from north pole to south pole of the outer ergosphere
        theta = np.linspace(0, np.pi, 100)

        #a turn around the outer ergosphere
        phi = np.linspace(0, 2 * np.pi, 100)

        #create dots
        theta, phi = np.meshgrid(theta, phi)

        #creating an array of radius based on each theta
        r = outer_ergosphere_radius(bh_mass, a_star, theta, False)

        #calling cartesian function to convert outer ergosphere coordinates from boyer-lindquist to cartesian
        x, y, z = cartesian(bh_mass, r, theta, phi, a_star)

        #creatig the black hole's outer ergosphere
        ax.plot_surface(x, y, z, color="deepskyblue", edgecolor="#2572A1", alpha=0.1, linewidth=0.1)

        #setting base for the coordinates of the polar informations
        polar_limit = outer_ergosphere_radius(bh_mass, a_star, radians(0), False)

        #setting base for the coordinates of the equatorial informations
        equatorial_limit = cartesian_radius(bh_mass, outer_ergosphere_radius(bh_mass, a_star, radians(90), False), a_star, np.radians(90))

        angle = np.radians(200)

        #setting coordinates for text and arrow for inner ergosphere cartesian equatorial radius
        y = polar_limit * np.cos(angle)
        z = polar_limit

        #drawing an arrow for outer ergosphere cartesian polar radius
        ax.plot([0, 0], [0, y], [z, z * 1.7], color="deepskyblue")

        #adding the text for outer ergosphere cartesian polar radius
        ax.text(0, y, z * 1.7, f"Outer Ergosphere Cartesian Polar Radius: " f"{length_unit_manager(polar_limit)}", path_effects=[pe.withStroke(linewidth=2, foreground="white")], color="deepskyblue", fontweight="bold", fontsize=10)


        angle = np.radians(180)

        #setting coordinates for text and arrow for inner ergosphere cartesian equatorial radius
        x = equatorial_limit * np.sin(angle)
        y = equatorial_limit * np.cos(angle)
        z = 0

        #drawing an arrow for outer ergosphere cartesian equatorial radius
        ax.plot([x, x * 1.5], [y, y * 1.5], [z, z], color="deepskyblue")

        #adding the text for outer ergosphere cartesian equatorial radius
        ax.text(x, y * 1.5, z, f"Outer Ergosphere Cartesian Equatorial Radius: " f"{length_unit_manager(equatorial_limit)}", path_effects=[pe.withStroke(linewidth=2, foreground="white")], color="deepskyblue", fontweight="bold", fontsize=10)

        ax.set_xlim(-equatorial_limit * 1.2, equatorial_limit * 1.2)
        ax.set_ylim(-equatorial_limit * 1.2, equatorial_limit * 1.2)
        ax.set_zlim(-equatorial_limit * 1.2, equatorial_limit * 1.2)

        ax.set_box_aspect([1, 1, 1])

        #show the drawings
        plt.pause(0.5)

#space-time structures
def space_time_draw(dimension, event_horizon_radius):
    #setting limits 
    limit = event_horizon_radius * 6

    flat_lines = []
    surfaces = []

    match dimension:
        case "2D":
            #setting limits
            x = np.linspace(-limit, limit, 40)
            y = np.linspace(-limit, limit, 40)
            #meshing dots
            x, y = np.meshgrid(x, y)

            z = np.zeros_like(x) - event_horizon_radius * 1.6

            #draw x lines
            for row in range(x.shape[0]):
                line, = ax.plot(x[row, :], y[row, :], z[row, :], color="white", linewidth=2, alpha=0.8)

                #saving x lines in the flat_lines
                flat_lines.append(line)

            #draw y lines
            for col in range(x.shape[1]):
                line, = ax.plot(x[:, col], y[:, col], z[:, col], color="white", linewidth=2, alpha=0.8)

                #saving y lines in the flat_lines
                flat_lines.append(line)

            #drawing and saving the surface
            surfaces = ax.plot_surface(x, y, z, color="#192841", linewidth=0, alpha=0.5)
        case "3D":
            line_alpha = 0.1
            line_color = "#505050"

            surface_alpha = 0.01
            surface_facecolor = "#586B8B"

            #setting limits
            x = np.linspace(-limit, limit, 20)
            y = np.linspace(-limit, limit, 20)
            z = np.linspace(-limit, limit, 20)

            #meshing dots
            x, y, z = np.meshgrid(x, y, z)

            #x
            x_lines = []
            x_surfaces = []

            #x lines
            for row_i in range(x.shape[1]):
                for row_j in range(x.shape[2]):
                    line = list(zip(x[:, row_i, row_j], y[:, row_i, row_j], z[:, row_i, row_j]))

                    #saving lines in the x_lines
                    x_lines.append(line)

                #x surfaces
                surface = [
                    (
                        x[row_i, 0, 0],
                        y[row_i, 0, 0],
                        z[row_i, 0, 0]
                    ),
                    (
                        x[row_i, -1, 0],
                        y[row_i, -1, 0],
                        z[row_i, -1, 0]
                    ),
                    (
                        x[row_i, -1, -1],
                        y[row_i, -1, -1],
                        z[row_i, -1, -1]
                    ),
                    (
                        x[row_i, 0, -1],
                        y[row_i, 0, -1],
                        z[row_i, 0, -1]
                    )
                ]

                #saving surfaces in the x_surfaces
                x_surfaces.append(surface)

            #all x lines becomes as one collection
            x_collection = Line3DCollection(x_lines, color=line_color, linewidth=1, alpha=line_alpha)

            #plot x lines
            ax.add_collection3d(x_collection)

            #saving x_collection in the flat_lines
            flat_lines.append(x_collection)

                
            #all x surfaces becomes as one collection
            x_surface_collection = Poly3DCollection(x_surfaces, facecolor=surface_facecolor, alpha=surface_alpha)

            #plot x surfaces
            ax.add_collection3d(x_surface_collection)

            #saving x_surface_collection in the surfaces
            surfaces.append(x_surface_collection)




            #y
            y_lines = []
            y_surfaces = []

            #y_lines
            for col_i in range(y.shape[1]):
                for col_j in range(y.shape[2]):
                    line = list(zip(x[col_i, :, col_j], y[col_i, :, col_j], z[col_i, :, col_j]))

                    #saving lines in the y_lines
                    y_lines.append(line)

                #y surfaces
                surface = [
                    (
                        x[0, col_i, 0],
                        y[0, col_i, 0],
                        z[0, col_i, 0]
                    ),
                    (
                        x[-1, col_i, 0],
                        y[-1, col_i, 0],
                        z[-1, col_i, 0]
                    ),
                    (
                        x[-1, col_i, -1],
                        y[-1, col_i, -1],
                        z[-1, col_i, -1]
                    ),
                    (
                        x[0, col_i, -1],
                        y[0, col_i, -1],
                        z[0, col_i, -1]
                    )
                ]
                
                #saving surfaces in the y_surfaces
                y_surfaces.append(surface)

            #all y lines becomes as one collection
            y_collection = Line3DCollection(y_lines, color=line_color, linewidth=1, alpha=line_alpha)

            #plot y lines
            ax.add_collection3d(y_collection)

            #saving x_collection in the flat_lines
            flat_lines.append(y_collection)

                
            #all y surfaces becomes as one collection
            y_surface_collection = Poly3DCollection(y_surfaces, facecolor=surface_facecolor, alpha=surface_alpha)

            #plot y surfaces
            ax.add_collection3d(y_surface_collection)

            #saving y_surface_collection in the surfaces
            surfaces.append(y_surface_collection)




            #z
            z_lines = []
            z_surfaces = []

            #z_lines
            for dep_i in range(z.shape[1]):
                for dep_j in range(z.shape[2]):
                    line = list(zip(x[dep_i, dep_j, :], y[dep_i, dep_j, :], z[dep_i, dep_j, :]))

                    #saving lines in the z_lines
                    z_lines.append(line)

                #z_surfaces
                surface = [
                    (
                        x[0, 0, dep_i],
                        y[0, 0, dep_i],
                        z[0, 0, dep_i]
                    ),
                    (
                        x[0, -1, dep_i],
                        y[0, -1, dep_i],
                        z[0, -1, dep_i]
                    ),
                    (
                        x[-1, -1, dep_i],
                        y[-1, -1, dep_i],
                        z[-1, -1, dep_i]
                    ),
                    (
                        x[-1, 0, dep_i],
                        y[-1, 0, dep_i],
                        z[-1, 0, dep_i]
                    )
                ]

                #saving surfaces in the z_surfaces
                z_surfaces.append(surface)

            #all z lines becomes as one collection
            z_collection = Line3DCollection(z_lines, color=line_color, linewidth=1, alpha=line_alpha)

            #plot z lines
            ax.add_collection3d(z_collection)

            #saving z_collection in the flat_lines
            flat_lines.append(z_collection)


            #all z surfaces becomes as one collection
            z_surface_collection = Poly3DCollection(z_surfaces, facecolor=surface_facecolor, alpha=surface_alpha)

            #plot z surfaces
            ax.add_collection3d(z_surface_collection)

            #saving z_surface_collection in the surfaces
            surfaces.append(z_surface_collection)




            # #top
            # surface = ax.plot_surface(x[:, :, 0], y[:, :, 0], z[:, :, 0], color="#192841", linewidth=0, alpha=0.5)
            # surfaces.append(surface)

            # #bottom
            # surface = ax.plot_surface(x[:, :, -1], y[:, :, -1], z[:, :, -1], color="#192841", linewidth=0, alpha=0.5)
            # surfaces.append(surface)

            # #front
            # surface = ax.plot_surface(x[:, 0, :], y[:, 0, :], z[:, 0, :], color="#192841", linewidth=0, alpha=0.5)
            # surfaces.append(surface)

            # #back
            # surface = ax.plot_surface(x[:, -1, :], y[:, -1, :], z[:, -1, :], color="#192841", linewidth=0, alpha=0.5)
            # surfaces.append(surface)

            # #right
            # surface = ax.plot_surface(x[0, :, :], y[0, :, :], z[0, :, :], color="#192841", linewidth=0, alpha=0.5)
            # surfaces.append(surface)

            # #left
            # surface = ax.plot_surface(x[-1, :, :], y[-1, :, :], z[-1, :, :], color="#192841", linewidth=0, alpha=0.5)
            # surfaces.append(surface)

    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    ax.set_zlim(-limit, limit)
    
    ax.set_box_aspect([1, 1, 1])

    #show the drawings
    plt.pause(0.5)

    return x, y, z, flat_lines, surfaces

def bend_space_time(x ,y, z, flat_lines, surfaces, dimension, event_horizon_radius):
    #removing flat_lines
    for line in flat_lines:
        line.remove()

    match dimension:
        case "2D":
            surfaces.remove()

            flat_lines = []
            surfaces = []

            #calculating the distance of each dot from the (0, 0)
            r = np.sqrt(x ** 2 + y ** 2)

            #main formula for setting z depends to the horizon's radius
            epsilon = event_horizon_radius * 0.001

            z = (event_horizon_radius * np.log(np.maximum(r - event_horizon_radius, epsilon) / event_horizon_radius)) - event_horizon_radius * 1.6

            #draw horizontal lines
            for row in range(x.shape[0]):
                line, = ax.plot(x[row, :], y[row, :], z[row, :], color="white", linewidth=2, alpha=0.8)

                #saving horizontal lines in the flat_lines
                flat_lines.append(line)

            #draw vertical lines
            for col in range(x.shape[1]):
                line, = ax.plot(x[:, col], y[:, col], z[:, col], color="white", linewidth=2, alpha=0.8)

                #saving vertical lines in the flat_lines
                flat_lines.append(line)    

            #drawing and saving the surface
            surface = ax.plot_surface(x, y, z, color="#192841", linewidth=0, alpha=0.2, shade=False)  
        case "3D":
            line_alpha = 0.1
            surface_alpha = 0.01

            #removing surface
            for surface in surfaces:
                surface.remove()

            flat_lines = []
            surfaces = []

            #calculating the distance of each dot from the (0, 0, 0)
            r = np.sqrt(x ** 2 + y ** 2 + z ** 2)

            #checking if r equals zero (preventing zero division)
            r_safe = np.where(r > 0, r, 1)

            #main formula for calculating the bend factor
            epsilon = event_horizon_radius * 0.001

            bend_factor = (event_horizon_radius * np.log(np.maximum(r - event_horizon_radius, epsilon) / event_horizon_radius)) - event_horizon_radius * 1.6

            #calculating the new radial distance
            new_r = np.maximum(r + bend_factor, epsilon)

            #calculating the radial scale
            scale = new_r / r_safe

            #applying the bend factor to coordinates
            x = x * scale
            y = y * scale
            z = z * scale

            bend_color_strength = np.clip(np.abs(bend_factor) / (event_horizon_radius * 5), 0, 1) * 255

            #creating r1 and r2 empty array like 3D-shape of bend_color_strength
            r1 = np.empty(bend_color_strength.shape, dtype=str)
            r2 = np.empty(bend_color_strength.shape, dtype=str)
            line_colors = np.empty(bend_color_strength.shape, dtype=object)
            surface_facecolors = np.empty(bend_color_strength.shape, dtype=object)

            red_strength_quotient = bend_color_strength // 16

            red_strength_remainder = bend_color_strength % 16

            for i in range(bend_color_strength.shape[0]):
                for j in range(bend_color_strength.shape[1]):
                    for k in range(bend_color_strength.shape[2]):
                        r1[i][j][k] = color_hex_list[int(red_strength_quotient[i][j][k])][1]
                        r2[i][j][k] = color_hex_list[int(red_strength_remainder[i][j][k])][1]

                        line_colors[i][j][k] = "#" + (r1)[i][j][k] + (r2)[i][j][k] + "1010"
                        surface_facecolors[i][j][k] = "#" + (r1)[i][j][k] + (r2)[i][j][k] + "1010"
            
            #x
            x_lines = []
            x_lines_colors = []
            x_surfaces = []
            x_surfaces_colors = []

            #x lines
            for row_i in range(x.shape[1]):
                for row_j in range(x.shape[2]):
                    line = list(zip(x[:, row_i, row_j], y[:, row_i, row_j], z[:, row_i, row_j]))

                    #saving lines in the x_lines
                    x_lines.append(line)

                    #finding the closest point to the center and extracting it's color
                    x_lines_colors.append(line_colors[np.argmin(r[:, row_i, row_j])][row_i][row_j])

                #x surfaces
                surface = [
                    (
                        x[row_i, 0, 0],
                        y[row_i, 0, 0],
                        z[row_i, 0, 0]
                    ),
                    (
                        x[row_i, -1, 0],
                        y[row_i, -1, 0],
                        z[row_i, -1, 0]
                    ),
                    (
                        x[row_i, -1, -1],
                        y[row_i, -1, -1],
                        z[row_i, -1, -1]
                    ),
                    (
                        x[row_i, 0, -1],
                        y[row_i, 0, -1],
                        z[row_i, 0, -1]
                    )
                ]

                #saving surfaces in the x_surfaces
                x_surfaces.append(surface)

                #finding the closest point to the center and extracting it's color
                x_surfaces_colors.append(surface_facecolors[row_i, x.shape[1] // 2, x.shape[2] // 2])

            #all x lines becomes as one collection
            x_collection = Line3DCollection(x_lines, color=x_lines_colors, linewidth=0.5, alpha=line_alpha)

            #plot x lines
            ax.add_collection3d(x_collection)

            #saving x_collection in the flat_lines
            flat_lines.append(x_collection)

                
            #all x surfaces becomes as one collection
            x_surface_collection = Poly3DCollection(x_surfaces, facecolor=x_surfaces_colors, alpha=surface_alpha)

            #plot x surfaces
            ax.add_collection3d(x_surface_collection)

            #saving x_surface_collection in the surfaces
            surfaces.append(x_surface_collection)




            #y
            y_lines = []
            y_lines_colors = []
            y_surfaces = []
            y_surfaces_colors = []

            #y_lines
            for col_i in range(y.shape[1]):
                for col_j in range(y.shape[2]):
                    line = list(zip(x[col_i, :, col_j], y[col_i, :, col_j], z[col_i, :, col_j]))

                    #saving lines in the y_lines
                    y_lines.append(line)

                    #finding the closest point to the center and extracting it's color
                    y_lines_colors.append(line_colors[col_i][np.argmin(r[col_i, :, col_j])][col_j])

                #y surfaces
                surface = [
                    (
                        x[0, col_i, 0],
                        y[0, col_i, 0],
                        z[0, col_i, 0]
                    ),
                    (
                        x[-1, col_i, 0],
                        y[-1, col_i, 0],
                        z[-1, col_i, 0]
                    ),
                    (
                        x[-1, col_i, -1],
                        y[-1, col_i, -1],
                        z[-1, col_i, -1]
                    ),
                    (
                        x[0, col_i, -1],
                        y[0, col_i, -1],
                        z[0, col_i, -1]
                    )
                ]
                
                #saving surfaces in the y_surfaces
                y_surfaces.append(surface)

                #finding the closest point to the center and extracting it's color
                y_surfaces_colors.append(surface_facecolors[y.shape[0] // 2, col_i, y.shape[2] // 2])

            #all y lines becomes as one collection
            y_collection = Line3DCollection(y_lines, color=y_lines_colors, linewidth=0.5, alpha=line_alpha)

            #plot y lines
            ax.add_collection3d(y_collection)

            #saving x_collection in the flat_lines
            flat_lines.append(y_collection)

                
            #all y surfaces becomes as one collection
            y_surface_collection = Poly3DCollection(y_surfaces, facecolor=y_surfaces_colors, alpha=surface_alpha)

            #plot y surfaces
            ax.add_collection3d(y_surface_collection)

            #saving y_surface_collection in the surfaces
            surfaces.append(y_surface_collection)




            #z
            z_lines = []
            z_lines_colors = []
            z_surfaces = []
            z_surfaces_colors = []

            #z_lines
            for dep_i in range(z.shape[1]):
                for dep_j in range(z.shape[2]):
                    line = list(zip(x[dep_i, dep_j, :], y[dep_i, dep_j, :], z[dep_i, dep_j, :]))

                    #saving lines in the z_lines
                    z_lines.append(line)

                    z_lines_colors.append(line_colors[dep_i][dep_j][np.argmin(r[dep_i, dep_j, :])])

                #z_surfaces
                surface = [
                    (
                        x[0, 0, dep_i],
                        y[0, 0, dep_i],
                        z[0, 0, dep_i]
                    ),
                    (
                        x[0, -1, dep_i],
                        y[0, -1, dep_i],
                        z[0, -1, dep_i]
                    ),
                    (
                        x[-1, -1, dep_i],
                        y[-1, -1, dep_i],
                        z[-1, -1, dep_i]
                    ),
                    (
                        x[-1, 0, dep_i],
                        y[-1, 0, dep_i],
                        z[-1, 0, dep_i]
                    )
                ]

                #saving surfaces in the z_surfaces
                z_surfaces.append(surface)

                #finding the closest point to the center and extracting it's color
                z_surfaces_colors.append(surface_facecolors[z.shape[0] // 2, z.shape[1] // 2, dep_i])

            #all z lines becomes as one collection
            z_collection = Line3DCollection(z_lines, color=z_lines_colors, linewidth=0.5, alpha=line_alpha)

            #plot z lines
            ax.add_collection3d(z_collection)

            #saving z_collection in the flat_lines
            flat_lines.append(z_collection)


            #all z surfaces becomes as one collection
            z_surface_collection = Poly3DCollection(z_surfaces, facecolor=z_surfaces_colors, alpha=surface_alpha)

            #plot z surfaces
            ax.add_collection3d(z_surface_collection)

            #saving z_surface_collection in the surfaces
            surfaces.append(z_surface_collection)
                        
    #show the drawings
    plt.pause(0.5)

    return x, y, z, flat_lines, surfaces

def twist_space_time(bh_mass, a_star, x, y, z, flat_lines, surfaces, dimension, event_horizon_radius):
    if a_star != 0:
        #removing flat_lines
        for line in flat_lines:
            line.remove()
        match dimension:
            case "2D":
                #removing surface
                for surface in surfaces:
                    surface.remove()

                flat_lines = []
                surface = []

                #creating an array of distances of every point in the space-time network
                r = np.sqrt((x ** 2) + (y ** 2))

                #calculating frame dragging angular velocity for each point and saving them into an array of omegas
                omega = frame_dragging_angular_velocity(bh_mass, a_star, r, gravitational_radius(bh_mass), np.pi / 2, False)

                #calculating twist strength for each omega
                twist_strength = omega * event_horizon_radius / c * 5

                x_twist = x * np.cos(twist_strength) - y * np.sin(twist_strength)
                y_twist = x * np.sin(twist_strength) + y * np.cos(twist_strength)
                z_twist = z

                #draw horizontal lines
                for row in range(x.shape[0]):
                    line, = ax.plot(x_twist[row, :], y_twist[row, :], z_twist[row, :], color="black", linewidth=1, alpha=0.2)

                    #saving horizontal lines in the flat_lines
                    flat_lines.append(line)

                #draw vertical lines
                for col in range(x.shape[1]):
                    line, = ax.plot(x_twist[:, col], y_twist[:, col], z_twist[:, col], color="black", linewidth=1, alpha=0.2)

                    #saving vertical lines in the flat_lines
                    flat_lines.append(line)    

                #drawing and saving the surface
                surface = ax.plot_surface(x_twist, y_twist, z_twist, color="#192841", linewidth=0, alpha=0.2)   
            case "3D":
                line_alpha = 0.1
                surface_alpha = 0.01

                #removing surface
                for surface in surfaces:
                    surface.remove()

                flat_lines = []
                surfaces = []

                #creating an array of distances of every point in the space-time network
                r = np.sqrt((x ** 2) + (y ** 2) + (z ** 2))

                #calculating frame dragging angular velocity for each point and saving them into an array of omegas
                omega = frame_dragging_angular_velocity(bh_mass, a_star, r, gravitational_radius(bh_mass), np.pi / 2, False)

                #calculating twist strength for each omega
                twist_strength = omega * event_horizon_radius / c * 5

                x_twist = x * np.cos(twist_strength) - y * np.sin(twist_strength)
                y_twist = x * np.sin(twist_strength) + y * np.cos(twist_strength)
                z_twist = z

                twist_color_strength = np.clip(np.abs(twist_strength) / 5, 0, 1) * 255

                #creating r1 and r2 empty array like 3D-shape of bend_color_strength
                r1 = np.empty(twist_color_strength.shape, dtype=str)
                r2 = np.empty(twist_color_strength.shape, dtype=str)
                line_colors = np.empty(twist_color_strength.shape, dtype=object)
                surface_facecolors = np.empty(twist_color_strength.shape, dtype=object)

                red_strength_quotient = twist_color_strength // 16

                red_strength_remainder = twist_color_strength % 16

                for i in range(twist_color_strength.shape[0]):
                    for j in range(twist_color_strength.shape[1]):
                        for k in range(twist_color_strength.shape[2]):
                            r1[i][j][k] = color_hex_list[int(red_strength_quotient[i][j][k])][1]
                            r2[i][j][k] = color_hex_list[int(red_strength_remainder[i][j][k])][1]

                            line_colors[i][j][k] = "#" + (r1)[i][j][k] + (r2)[i][j][k] + "1010"
                            surface_facecolors[i][j][k] = "#243A" + (r1)[i][j][k] + (r2)[i][j][k]
                #x
                x_lines = []
                x_lines_colors = []
                x_surfaces = []
                x_surfaces_colors = []

                #x lines
                for row_i in range(x.shape[1]):
                    for row_j in range(x.shape[2]):
                        line = list(zip(x_twist[:, row_i, row_j], y_twist[:, row_i, row_j], z_twist[:, row_i, row_j]))

                        #saving lines in the x_lines
                        x_lines.append(line)

                        #finding the closest point to the center and extracting it's color
                        x_lines_colors.append(line_colors[np.argmin(r[:, row_i, row_j])][row_i][row_j])

                    #x surfaces
                    surface = [
                        (
                            x_twist[row_i, 0, 0],
                            y_twist[row_i, 0, 0],
                            z_twist[row_i, 0, 0]
                        ),
                        (
                            x_twist[row_i, -1, 0],
                            y_twist[row_i, -1, 0],
                            z_twist[row_i, -1, 0]
                        ),
                        (
                            x_twist[row_i, -1, -1],
                            y_twist[row_i, -1, -1],
                            z_twist[row_i, -1, -1]
                        ),
                        (
                            x_twist[row_i, 0, -1],
                            y_twist[row_i, 0, -1],
                            z_twist[row_i, 0, -1]
                        )
                    ]

                    #saving surfaces in the x_surfaces
                    x_surfaces.append(surface)

                    #finding the closest point to the center and extracting it's color
                    x_surfaces_colors.append(surface_facecolors[row_i, x.shape[1] // 2, x.shape[2] // 2])

                #all x lines becomes as one collection
                x_collection = Line3DCollection(x_lines, color=x_lines_colors, linewidth=0.5, alpha=line_alpha)

                #plot x lines
                ax.add_collection3d(x_collection)

                #saving x_collection in the flat_lines
                flat_lines.append(x_collection)

                    
                #all x surfaces becomes as one collection
                x_surface_collection = Poly3DCollection(x_surfaces, facecolor=x_surfaces_colors, alpha=surface_alpha)

                #plot x surfaces
                ax.add_collection3d(x_surface_collection)

                #saving x_surface_collection in the surfaces
                surfaces.append(x_surface_collection)




                #y
                y_lines = []
                y_lines_colors = []
                y_surfaces = []
                y_surfaces_colors = []

                #y_lines
                for col_i in range(y.shape[1]):
                    for col_j in range(y.shape[2]):
                        line = list(zip(x_twist[col_i, :, col_j], y_twist[col_i, :, col_j], z_twist[col_i, :, col_j]))

                        #saving lines in the y_lines
                        y_lines.append(line)

                        #finding the closest point to the center and extracting it's color
                        y_lines_colors.append(line_colors[col_i][np.argmin(r[col_i, :, col_j])][col_j])

                    #y surfaces
                    surface = [
                        (
                            x_twist[0, col_i, 0],
                            y_twist[0, col_i, 0],
                            z_twist[0, col_i, 0]
                        ),
                        (
                            x_twist[-1, col_i, 0],
                            y_twist[-1, col_i, 0],
                            z_twist[-1, col_i, 0]
                        ),
                        (
                            x_twist[-1, col_i, -1],
                            y_twist[-1, col_i, -1],
                            z_twist[-1, col_i, -1]
                        ),
                        (
                            x_twist[0, col_i, -1],
                            y_twist[0, col_i, -1],
                            z_twist[0, col_i, -1]
                        )
                    ]
                    
                    #saving surfaces in the y_surfaces
                    y_surfaces.append(surface)

                    #finding the closest point to the center and extracting it's color
                    y_surfaces_colors.append(surface_facecolors[y.shape[0] // 2, col_i, y.shape[2] // 2])

                #all y lines becomes as one collection
                y_collection = Line3DCollection(y_lines, color=y_lines_colors, linewidth=0.5, alpha=line_alpha)

                #plot y lines
                ax.add_collection3d(y_collection)

                #saving x_collection in the flat_lines
                flat_lines.append(y_collection)

                    
                #all y surfaces becomes as one collection
                y_surface_collection = Poly3DCollection(y_surfaces, facecolor=y_surfaces_colors, alpha=surface_alpha)

                #plot y surfaces
                ax.add_collection3d(y_surface_collection)

                #saving y_surface_collection in the surfaces
                surfaces.append(y_surface_collection)




                #z
                z_lines = []
                z_lines_colors = []
                z_surfaces = []
                z_surfaces_colors = []

                #z_lines
                for dep_i in range(z.shape[1]):
                    for dep_j in range(z.shape[2]):
                        line = list(zip(x_twist[dep_i, dep_j, :], y_twist[dep_i, dep_j, :], z_twist[dep_i, dep_j, :]))

                        #saving lines in the z_lines
                        z_lines.append(line)

                        z_lines_colors.append(line_colors[dep_i][dep_j][np.argmin(r[dep_i, dep_j, :])])

                    #z_surfaces
                    surface = [
                        (
                            x_twist[0, 0, dep_i],
                            y_twist[0, 0, dep_i],
                            z_twist[0, 0, dep_i]
                        ),
                        (
                            x_twist[0, -1, dep_i],
                            y_twist[0, -1, dep_i],
                            z_twist[0, -1, dep_i]
                        ),
                        (
                            x_twist[-1, -1, dep_i],
                            y_twist[-1, -1, dep_i],
                            z_twist[-1, -1, dep_i]
                        ),
                        (
                            x_twist[-1, 0, dep_i],
                            y_twist[-1, 0, dep_i],
                            z_twist[-1, 0, dep_i]
                        )
                    ]

                    #saving surfaces in the z_surfaces
                    z_surfaces.append(surface)

                    #finding the closest point to the center and extracting it's color
                    z_surfaces_colors.append(surface_facecolors[z.shape[0] // 2, z.shape[2] // 2, dep_i])

                #all z lines becomes as one collection
                z_collection = Line3DCollection(z_lines, color=z_lines_colors, linewidth=0.5, alpha=line_alpha)

                #plot z lines
                ax.add_collection3d(z_collection)

                #saving z_collection in the flat_lines
                flat_lines.append(z_collection)


                #all z surfaces becomes as one collection
                z_surface_collection = Poly3DCollection(z_surfaces, facecolor=z_surfaces_colors, alpha=surface_alpha)

                #plot z surfaces
                ax.add_collection3d(z_surface_collection)

                #saving z_surface_collection in the surfaces
                surfaces.append(z_surface_collection)

    #show the drawings
    plt.pause(0.5)

    return x_twist, y_twist, z_twist, flat_lines, surfaces

#orbital structures

#photon orbitals
def photon_sphere_radius(bh_type, bh_mass, a_star, direction, theta, show):
    if a_star == 0:
        phs_radius = 3 * gravitational_radius(bh_mass)
        if show == True:
            print(Fore.LIGHTGREEN_EX + f"Photon Sphere Radius at {round(np.degrees(theta))}°: {Fore.RESET} {length_unit_manager(phs_radius)}", end="\n")
    else:
        if a_star == 1 and direction == "prograde":
            a_star = 0.999999999999999

        m = gravitational_radius(bh_mass)

        roots = []

        #creating 10000 possible radiuses to check if it's the photon ring or not
        r_possible_values = np.linspace(np.nextafter(outer_event_horizon_radius(bh_type, bh_mass, a_star, theta, False), np.inf), m * 4, 10000)

        tolerance = 1e-10

        for r1, r2 in zip(r_possible_values[:-1], r_possible_values[1:]):

            f1 = photon_equation(bh_mass, a_star, r1, direction, theta)
            f2 = photon_equation(bh_mass, a_star, r2, direction, theta)
        
            if not np.isfinite(f1) or not np.isfinite(f2):
                continue

            if abs(f1) < tolerance:
                root = r1
            elif abs(f2) < tolerance:
                root = r2
            elif f1 * f2 < 0:
                root = opt.brentq(
                    lambda r: photon_equation(bh_mass, a_star, r, direction, theta),
                    r1,
                    r2
                )

            else:
                continue

            if not roots or abs(root - roots[-1]) > tolerance:
                roots.append(root)

        match direction:
            case "prograde":
                phs_radius = min(roots)
                if show == True:
                    print(Fore.LIGHTGREEN_EX + f"Co-Rotating Photon Sphere Boyer-Lindquist Radius at {round(np.degrees(theta))}°: {Fore.RESET} {length_unit_manager(phs_radius)}", end="\n")
                    print(Fore.LIGHTCYAN_EX + f"Co-Rotating Photon Sphere Cartesian Radius at {round(np.degrees(theta))}°: {Fore.RESET} {length_unit_manager(cartesian_radius(bh_mass, phs_radius, a_star, theta))}", end="\n")
            case "retrograde":
                phs_radius = max(roots)
                if show == True:
                    print(Fore.LIGHTGREEN_EX + f"Counter-Rotating Photon Sphere Boyer-Lindquist Radius at {round(np.degrees(theta))}°: {Fore.RESET} {length_unit_manager(phs_radius)}", end="\n")
                    print(Fore.LIGHTCYAN_EX + f"Counter-Rotating Photon Sphere Cartesian Radius at {round(np.degrees(theta))}°: {Fore.RESET} {length_unit_manager(cartesian_radius(bh_mass, phs_radius, a_star, theta))}", end="\n")

        if show == True:
            fda_velocity = frame_dragging_angular_velocity(bh_mass, a_star, phs_radius, gravitational_radius(bh_mass), theta, False)
            fdl_velocity = frame_dragging_linear_velocity(bh_mass, a_star, fda_velocity, phs_radius, theta)
            print(Fore.RED + f"Frame-dragging coordinate linear speed at {round(np.degrees(theta))}°: {Fore.RESET} {speed_unit_manager(fdl_velocity)}", end="\n")

    if show == True:
        print("         ********************            ") 
        print(Fore.LIGHTRED_EX + f"it's about: {Fore.RESET} {length_unit_manager(cartesian_radius(bh_mass, phs_radius, a_star, theta) - cartesian_radius(bh_mass, outer_event_horizon_radius(bh_type, bh_mass, a_star, theta, False), a_star, theta))} {Fore.LIGHTRED_EX} from the Outer Event Horizon.", end="\n")
        print(Fore.RED + "----------------------------------------" + Fore.RESET)
    return phs_radius

def photon_ring_draw(bh_type, bh_mass, a_star, direction, theta):
    #setting the amount of the ring inclination
    inclination = np.pi / 2 - theta

    #a turn around the ring
    phi = np.linspace(0, 2 * np.pi, 100)

    if bh_type == "schw" or bh_type == "1":
        #getting radius from photon sphere function
        r = photon_sphere_radius(bh_type, bh_mass, a_star, "", theta, False)

        #converting angles and radii into coordinates
        x, y, z = cartesian(bh_mass, r, np.pi / 2, phi, a_star)

        #rotate the ring around the x-axis
        y_rot = y * np.cos(inclination) - z * np.sin(inclination)
        z_rot = y * np.sin(inclination) + z * np.cos(inclination)

        y = y_rot
        z = z_rot

        text = "Photon Ring Radius at"

        #creatig the black hole's photon sphere
        ax.plot(x, y, z, color="red", linewidth=2)

        x_text = r * np.sin(theta)
        y_text = r * np.cos(theta)
        z_text = 0

        # rotate text position around x-axis
        y_rot = y_text * np.cos(inclination) - z_text * np.sin(inclination)
        z_rot = y_text * np.sin(inclination) + z_text * np.cos(inclination)

        y_text = y_rot
        z_text = z_rot

        #drawing the arrow for photon sphere
        ax.plot([x_text, x_text * 1.3], [y_text, y_text * 1.3], [z_text, z_text], color="red")

        #adding the text for photon sphere
        ax.text(x_text * 1.3, y_text * 1.3, z_text, f"{text} {round(np.degrees(theta))}°: " f"{length_unit_manager(r)}", color="red", path_effects=[pe.withStroke(linewidth=2, foreground="white")], fontsize=10)
    elif bh_type == "kerr" or bh_type == "2":
        match direction:
            case "prograde":
                text = "Co-Rotating Photon Ring Cartesian Radius at"
                color = "red"
            case "retrograde":
                text = "Counter-Rotating Photon Ring Cartesian Radius at"
                color = "#ff7f3e"

        #getting radius from photon sphere function
        r = photon_sphere_radius(bh_type, bh_mass, a_star, direction, theta, False)

        #creating the photon ring
        if np.isclose(theta, 0) or np.isclose(theta, np.pi):
            ring_radius = cartesian_radius(bh_mass, r, a_star, theta)
            x = ring_radius * np.cos(phi)
            y = np.zeros_like(phi)
            z = ring_radius * np.sin(phi)
        else:
            x ,y, z = cartesian(bh_mass, r, np.pi / 2, phi, a_star)

            #rotate the ring around the x-axis
            y_rot = y * np.cos(inclination) - z * np.sin(inclination)
            z_rot = y * np.sin(inclination) + z * np.cos(inclination)

            y = y_rot
            z = z_rot

        #creatig the black hole's photon ring
        ax.plot(x, y, z, color=color, linewidth=2)

        #converting angles and radii into coordinates for texts
        x_text, y_text, z_text = cartesian(bh_mass, r, np.pi / 2, theta * 2, a_star)

        # rotate text position around x-axis
        y_rot = y_text * np.cos(inclination) - z_text * np.sin(inclination)
        z_rot = y_text * np.sin(inclination) + z_text * np.cos(inclination)

        y_text = y_rot
        z_text = z_rot
        
        #drawing the arrow for photon ring
        ax.plot([x_text, x_text * 1.6], [y_text, y_text * 1.6], [z_text, z_text * 1.3], color=color)

        #adding space-time velocity at photon ring
        fda_velocity = frame_dragging_angular_velocity(bh_mass, a_star, r, gravitational_radius(bh_mass), theta, False)
        fdl_velocity = frame_dragging_linear_velocity(bh_mass, a_star, fda_velocity, r, theta)

        #adding the text for photon ring
        ax.text(x_text * 1.6, y_text * 1.6, z_text * 1.3, f"{text} {round(np.degrees(theta))}°: " f"{length_unit_manager(cartesian_radius(bh_mass, r, a_star, theta))} \n Frame-dragging coordinate linear speed: {speed_unit_manager(fdl_velocity)}", color=color, path_effects=[pe.withStroke(linewidth=2, foreground="white")], fontsize=10)

    #show the drawings
    plt.pause(1)

def photon_sphere_draw(bh_type, bh_mass, a_star, direction):

    #points from north pole to south pole
    theta = np.linspace(0, np.pi, 100)

    #a turn around the sphere
    phi = np.linspace(0, 2 * np.pi, 100)

    #creating a 2D array of theta and phi and storing them
    theta_grid, phi_grid = np.meshgrid(theta, phi)

    #for a schwarzschild black hole
    if bh_type in ["schw", "1"]:
        #setting the radius or distance to the center
        r_val = photon_sphere_radius(bh_type, bh_mass, a_star, direction, "", False)
        #setting the shape of the array of the distances like theta_grid or phi_grid
        r_grid = np.full_like(theta_grid, r_val)

        #converting arrays of angles and radii into coordinates for plot
        x, y, z = cartesian(bh_mass, r_grid, theta_grid, phi_grid, 0)
        #plot the surface 
        photon_surface = ax.plot_surface(x, y, z, color="red", alpha=0.1, edgecolor="darkred")
    elif bh_type in ["kerr", "2"]:
        #creating a list of radii or distances with the length of the theta and storing each radius based on each angle
        r_by_theta = []
        for t in theta:
            r = photon_sphere_radius(bh_type, bh_mass, a_star, direction, t, False)
            r_by_theta.append(r)

        #creating copies of the array of radii phi times in the rows and one time in the column
        r_grid = np.tile(r_by_theta, (len(phi), 1))

        #converting arrays of angles and radii into coordinates for plot
        x, y, z = cartesian(bh_mass, r_grid, theta_grid, phi_grid, a_star)
        
        if direction == "prograde":
            color = "red"
            edgecolor = "darkred"
        elif direction == "retrograde":
            color = "#ff7f3e"
            edgecolor = "#72391c"

        #plot the surface
        photon_surface = ax.plot_surface(x, y, z, color=color, alpha=0.1, edgecolor=edgecolor)

    plt.pause(0.5)
    return photon_surface

#isco orbitals
def isco_sphere_radius(bh_type, bh_mass, a_star, direction, theta, show):
    outer_ev = outer_event_horizon_radius(bh_type, bh_mass, a_star, theta, False)
    if bh_type in ["schw", "1"]:
        isco_radius = (gravitational_radius(bh_mass) * 6)
        if show == True:
            print(Fore.LIGHTGREEN_EX + f"Safe Orbital For Matter(ISCO) Radius at {round(np.degrees(theta))}°: {Fore.RESET} {length_unit_manager(isco_radius)}", end="\n")
            print(Fore.LIGHTRED_EX + f"it's about: {Fore.RESET} {length_unit_manager(isco_radius-outer_ev)} {Fore.LIGHTRED_EX} from the Outer Event Horizon.", end="\n")
            print(Fore.RED + "----------------------------------------" + Fore.RESET)

    elif bh_type in ["kerr", "2"]:
        if a_star == 1 and direction == "prograde":
            a_star = 0.999999999999999

        m = gravitational_radius(bh_mass)

        roots = []

        #creating 20000 possible radii to check if it's the ISCO ring or not
        r_possible_values = np.linspace(np.nextafter(outer_event_horizon_radius(bh_type, bh_mass, a_star, theta, False), np.inf), m * 9, 20000)

        tolerance = 1e-10

        for r1, r2 in zip(r_possible_values[:-1], r_possible_values[1:]):
            match theta:
                case _ if not np.isclose(theta, 0) and not np.isclose(theta, np.pi):
                    f1 = isco_equation(bh_mass, a_star, r1, direction, theta)
                    f2 = isco_equation(bh_mass, a_star, r2, direction, theta)
                case _ if np.isclose(theta, 0) or np.isclose(theta, np.pi):
                    f1 = polar_isco_equation(bh_mass, a_star, r1)
                    f2 = polar_isco_equation(bh_mass, a_star, r2)

            if not np.isfinite(f1) or not np.isfinite(f2):
                continue

            if abs(f1) < tolerance:
                root = r1
            elif abs(f2) < tolerance:
                root = r2
            elif f1 * f2 < 0:
                match theta:
                    case _ if not np.isclose(theta, 0) and not np.isclose(theta, np.pi):
                        root = opt.brentq(
                            lambda r: isco_equation(bh_mass, a_star, r, direction, theta),
                            r1,
                            r2
                        )
                    case _ if np.isclose(theta, 0) or np.isclose(theta, np.pi):
                        root = opt.brentq(
                            lambda r: polar_isco_equation(bh_mass, a_star, r),
                            r1,
                            r2
                        )
            else:
                continue

            if not roots or abs(root - roots[-1]) > tolerance:
                roots.append(root)

        match direction:
            case "prograde":
                isco_radius = min(roots)
                if show == True:
                    print(Fore.LIGHTGREEN_EX + f"Prograde ISCO Boyer-Lindquist Radius at {round(np.degrees(theta))}°: {Fore.RESET} {length_unit_manager(isco_radius)}", end="\n")
                    print(Fore.LIGHTCYAN_EX + f"Prograde ISCO Cartesian Radius at {round(np.degrees(theta))}°: {Fore.RESET} {length_unit_manager(cartesian_radius(bh_mass, isco_radius, a_star, theta))}", end="\n")
            case "retrograde":
                isco_radius = max(roots)
                if show == True:
                    print(Fore.LIGHTGREEN_EX + f"Retrograde ISCO Boyer-Lindquist Radius at {round(np.degrees(theta))}°: {Fore.RESET} {length_unit_manager(isco_radius)}", end="\n")
                    print(Fore.LIGHTCYAN_EX + f"Retrograde ISCO Cartesian Radius at {round(np.degrees(theta))}°: {Fore.RESET} {length_unit_manager(cartesian_radius(bh_mass, isco_radius, a_star, theta))}", end="\n")

        if show == True:
            fda_velocity = frame_dragging_angular_velocity(bh_mass, a_star, isco_radius, gravitational_radius(bh_mass), theta, False)
            fdl_velocity = frame_dragging_linear_velocity(bh_mass, a_star, fda_velocity, isco_radius, theta)
            print(Fore.RED + f"Frame-dragging coordinate linear speed at {round(np.degrees(theta))}°: {Fore.RESET} {speed_unit_manager(fdl_velocity)}", end="\n")

    if show == True:
        print("         ********************            ") 
        print(Fore.LIGHTRED_EX + f"you shoud orbit at least at a distance of {Fore.RESET} {length_unit_manager(isco_radius-outer_ev)} {Fore.LIGHTRED_EX} from the Outer Event Horizon.")
        print(Fore.RED + "----------------------------------------" + Fore.RESET)

    return isco_radius

def isco_ring_draw(bh_type, bh_mass, a_star, direction, theta):
    #setting the amount of the ring inclination
    inclination = np.pi / 2 - theta

    #a turn around the ISCO
    phi = np.linspace(0, 2 * np.pi, 100)
    
    #getting radius from ISCO function
    r = isco_sphere_radius(bh_type, bh_mass, a_star, direction, theta, False)

    if bh_type in ["schw", "1"]:
        #converting angles and radii into coordinates
        x, y, z = cartesian(bh_mass, r, np.pi / 2, phi, a_star)

        #rotate the ring around the x-axis
        y_rot = y * np.cos(inclination) - z * np.sin(inclination)
        z_rot = y * np.sin(inclination) + z * np.cos(inclination)

        y = y_rot
        z = z_rot

        text = "The ISCO Radius at"

        #creatig the black hole's ISCO shpere
        ax.plot(x, y, z, color="#4cbb17", linewidth=2)

        x_text = r * np.sin(theta)
        y_text = r * np.cos(theta)
        z_text = 0

        #rotate text position around x-axis
        y_rot = y_text * np.cos(inclination) - z_text * np.sin(inclination)
        z_rot = y_text * np.sin(inclination) + z_text * np.cos(inclination)

        y_text = y_rot
        z_text = z_rot

        #drawing the arrow for ISCO sphere
        ax.plot([x_text, x_text * 1.3], [y_text, y_text * 1.3], [z_text, z_text], color="#4cbb17")

        #adding the text for ISCO sphere
        ax.text(x_text * 1.3, y_text * 1.3, z_text, f"{text} {round(np.degrees(theta))}°: " f"{length_unit_manager(r)}", color="#4cbb17", path_effects=[pe.withStroke(linewidth=2, foreground="white")], fontsize=10)
    elif bh_type in ["kerr", "2"]:
        match direction:
            case "prograde":
                text = "Co-Rotating ISCO Radius at"
                color = "#4cbb17"
            case "retrograde":
                text = "Counter-Rotating ISCO Radius at"
                color = "#0b6623"

        #creating the ISCO ring
        if np.isclose(theta, 0) or np.isclose(theta, np.pi):
            ring_radius = cartesian_radius(bh_mass, r, a_star, theta)
            x = ring_radius * np.cos(phi)
            y = np.zeros_like(phi)
            z = ring_radius * np.sin(phi)
        else:
            #converting angles and radii into coordinates
            x, y, z = cartesian(bh_mass, r, np.pi / 2, phi, a_star)

            #rotate the ring around the x-axis
            y_rot = y * np.cos(inclination) - z * np.sin(inclination)
            z_rot = y * np.sin(inclination) + z * np.cos(inclination)

            y = y_rot
            z = z_rot

        #creatig the black hole's ISCO radius
        ax.plot(x, y, z, color=color, linewidth=2)

        #converting angles and radii into coordinates for texts
        x_text, y_text, z_text = cartesian(bh_mass, r, np.pi / 2, theta * 2, a_star)

        x_text = r * np.sin(theta)
        y_text = r * np.cos(theta)
        z_text = 0

        #rotate text position around x-axis
        y_rot = y_text * np.cos(inclination) - z_text * np.sin(inclination)
        z_rot = y_text * np.sin(inclination) + z_text * np.cos(inclination)

        y_text = y_rot
        z_text = z_rot

        #drawing the arrow for ISCO sphere
        ax.plot([x_text, x_text * 1.3], [y_text, y_text * 1.3], [z_text, z_text], color=color)
        
        #adding space-time velocity at ISCO
        fda_velocity = frame_dragging_angular_velocity(bh_mass, a_star, r, gravitational_radius(bh_mass), theta, False)
        fdl_velocity = frame_dragging_linear_velocity(bh_mass, a_star, fda_velocity, r, theta)

        #adding the text for co-rotating ISCO sphere
        ax.text(x_text * 1.2, y_text * 1.2, z_text * 1.1, f"{text} {round(np.degrees(theta))}°: " f"{length_unit_manager(r)} \n Frame Dragging Velocity: {speed_unit_manager(fdl_velocity)}", color=color, path_effects=[pe.withStroke(linewidth=2, foreground="white")], fontsize=10)

    #show the drawings
    plt.pause(1)

def isco_sphere_draw(bh_type, bh_mass, a_star, direction):
    #points from north pole to south pole
    theta = np.linspace(0, np.pi, 100)

    #a turn around the sphere
    phi = np.linspace(0, 2 * np.pi, 100)

    #creating a 2D array of theta and phi and storing them
    theta_grid, phi_grid = np.meshgrid(theta, phi)

    #for a schwarzschild black hole
    if bh_type in ["schw", "1"]:
        #setting the radius or distance to the center
        r_val = isco_sphere_radius(bh_type, bh_mass, a_star, direction, "", False)
        #setting the shape of the array of the distances like theta_grid or phi_grid
        r_grid = np.full_like(theta_grid, r_val)

        #converting arrays of angles and radii into coordinates for plot
        x, y, z = cartesian(bh_mass, r_grid, theta_grid, phi_grid, 0)
        #plot the surface 
        isco_surface = ax.plot_surface(x, y, z, color="#4cbb17", alpha=0.1, edgecolor="#2A670D")
    elif bh_type in ["kerr", "2"]:
        #creating a list of radii or distances with the length of the theta and storing each radius based on each angle
        r_by_theta = []
        for t in theta:
            r = isco_sphere_radius(bh_type, bh_mass, a_star, direction, t, False)
            r_by_theta.append(r)

        #creating copies of the array of radii phi times in the rows and one time in the column
        r_grid = np.tile(r_by_theta, (len(phi), 1))

        #converting arrays of angles and radii into coordinates for plot
        x, y, z = cartesian(bh_mass, r_grid, theta_grid, phi_grid, a_star)
        
        if direction == "prograde":
            color = "#4cbb17"
            edgecolor = "#2A670D"
        elif direction == "retrograde":
            color = "#0b6623"
            edgecolor = "#104d20"

        #plot the surface
        isco_surface = ax.plot_surface(x, y, z, color=color, alpha=0.1, edgecolor=edgecolor)

    plt.pause(0.5)
    return isco_surface

#object experience
def stationary_object_xp_info_print(bh_type, bh_mass, a_star, distance, s_polar_angle, s_azimuthal_angle):
    x_center, y_center, z_center = cartesian(bh_mass, distance, s_polar_angle, s_azimuthal_angle, a_star)

    outer_ev = outer_event_horizon_radius(bh_type, bh_mass, a_star, s_polar_angle, False)

    print(Fore.LIGHTMAGENTA_EX + f"Outer Event Horizon:", end="\n")
    print(Fore.LIGHTGREEN_EX + f"Outer Event Horizon Boyer-Lindquist Radius at the Polar Angle of {Fore.LIGHTRED_EX}{np.round(np.degrees(s_polar_angle))}°{Fore.LIGHTGREEN_EX}: {Fore.RESET}{length_unit_manager(outer_ev)}.", end="\n")
    print(Fore.LIGHTCYAN_EX + f"Outer Event Horizon Cartesian Radius at the Polar Angle of {Fore.LIGHTRED_EX}{np.round(np.degrees(s_polar_angle))}°{Fore.LIGHTGREEN_EX}: {Fore.RESET}{length_unit_manager(cartesian_radius(bh_mass, outer_ev, a_star, s_polar_angle))}{Fore.LIGHTGREEN_EX}.", end="\n")
    print(Fore.RED + "----------------------------------------" + Fore.RESET)

    outer_er = outer_ergosphere_radius(bh_mass, a_star, s_polar_angle, False)

    print(Fore.LIGHTMAGENTA_EX + f"Outer Ergosphere:", end="\n")
    print(Fore.LIGHTGREEN_EX + f"Outer Ergosphere Boyer-Lindquist Radiusat at the Polar Angle of {Fore.LIGHTRED_EX}{np.round(np.degrees(s_polar_angle))}°{Fore.LIGHTGREEN_EX}: {Fore.RESET}{length_unit_manager(outer_er)}.", end="\n")
    print(Fore.LIGHTCYAN_EX + f"Outer Ergosphere Cartesian Radiusat the Polar Angle of {Fore.LIGHTRED_EX}{np.round(np.degrees(s_polar_angle))}°{Fore.LIGHTGREEN_EX}: {Fore.RESET}{length_unit_manager(cartesian_radius(bh_mass, outer_er, a_star, s_polar_angle))}", end="\n")
    print(Fore.RED + "----------------------------------------" + Fore.RESET)

    object_cartesian_distance_from_bh_center = np.sqrt(x_center**2 + y_center**2 + z_center**2)

    print(Fore.LIGHTCYAN_EX + f"Object's Cartesian Distance from the Black Hole's Center: {Fore.RESET}{length_unit_manager(object_cartesian_distance_from_bh_center)}.", end="\n")

    object_cartesian_distance_from_horizon = object_cartesian_distance_from_bh_center - cartesian_radius(bh_mass, outer_ev, a_star, s_polar_angle)

    print(Fore.LIGHTCYAN_EX + f"Object's Cartesian Distance from the Outer Event Horizon: {Fore.RESET}{length_unit_manager(object_cartesian_distance_from_horizon)}.", end="\n")
    print(Fore.RED + "----------------------------------------" + Fore.RESET)

    object_cartesian_distance_from_ergosphere = object_cartesian_distance_from_bh_center - cartesian_radius(bh_mass, outer_er, a_star, s_polar_angle)

    if object_cartesian_distance_from_ergosphere >= 0:
        print(Fore.LIGHTCYAN_EX + f"Object's Cartesian Distance from the Outer Ergosphere:{Fore.RESET}{length_unit_manager(object_cartesian_distance_from_ergosphere)}.", end="\n")
    else:
        print(Fore.LIGHTCYAN_EX + f"Object's Cartesian Distance from the Outer Ergosphere:{Fore.RED} NEGATIVE!{Fore.RESET}{length_unit_manager(np.abs(object_cartesian_distance_from_ergosphere))}.", end="\n")
    print(Fore.RED + "----------------------------------------" + Fore.RESET)

    if bh_type in ["kerr", "2"] and object_cartesian_distance_from_ergosphere <= 0:
        print(f"{Fore.RED}Object CAN'T be Stationary inside the Black Hole's Ergosphere!!\n"
              f"{Fore.LIGHTGREEN_EX}Stationary Object at a Distance of {Fore.RESET}{length_unit_manager(object_cartesian_distance_from_bh_center)}{Fore.LIGHTGREEN_EX} from the Black Hole's Center.\n"
                f"and {Fore.RESET}{length_unit_manager(object_cartesian_distance_from_horizon)}{Fore.LIGHTGREEN_EX} from the Outer Horizon.")
    else:
        print(f"{Fore.LIGHTGREEN_EX}Stationary Object Can be at a Distance of {Fore.RESET}{length_unit_manager(object_cartesian_distance_from_bh_center)}{Fore.LIGHTGREEN_EX} from the Black Hole's Center.\n"
              f"and {Fore.RESET}{length_unit_manager(object_cartesian_distance_from_horizon)}{Fore.LIGHTGREEN_EX} from the Outer Horizon.")
    print(Fore.RED + "----------------------------------------" + Fore.RESET)

    if bh_type in ["schw", "1"]:
        print(f"{Fore.LIGHTGREEN_EX}No Frame Dragging for a Schwarzschild Black Hole.{Fore.RESET}")
    elif bh_type in ["kerr", "2"]:
        fda_velocity = frame_dragging_angular_velocity(bh_mass, a_star, distance, gravitational_radius(bh_mass), s_polar_angle, False)
        fdl_velocity = frame_dragging_linear_velocity(bh_mass, a_star, fda_velocity, distance, s_polar_angle)
        print(f"{Fore.LIGHTGREEN_EX}Frame-dragging coordinate linear speed: {Fore.RESET}{speed_unit_manager(fdl_velocity)}.")

    return (
        x_center, y_center, z_center,
        outer_ev,
        object_cartesian_distance_from_bh_center,
        object_cartesian_distance_from_horizon,
        object_cartesian_distance_from_ergosphere,
        fdl_velocity
    )

def stationary_object_xp_draw(bh_type, bh_mass, a_star, distance, s_polar_angle, s_azimuthal_angle):
    object_data = stationary_object_xp_info_print(bh_type, bh_mass, a_star, distance, s_polar_angle, s_azimuthal_angle)

    x_center, y_center, z_center = object_data[0], object_data[1], object_data[2]

    outer_ev = object_data[3]

    object_cartesian_distance_from_bh_center = object_data[4]

    object_cartesian_distance_from_horizon = object_data[5]

    object_cartesian_distance_from_ergosphere = object_data[6]

    if bh_type in ["kerr", "2"] and object_cartesian_distance_from_ergosphere <= 0:
        edgecolor = "red"
        text = f"Object CAN'T be Stationary inside the Black Hole's Ergosphere!!\n Stationary Object at a Distance of {length_unit_manager(object_cartesian_distance_from_bh_center)} from the Black Hole's Center.\n and {length_unit_manager(object_cartesian_distance_from_horizon)} from the Outer Horizon."
    else:
        edgecolor = "green"
        text = f"Stationary Object at a Distance of {length_unit_manager(object_cartesian_distance_from_bh_center)} from the Black Hole's Center.\n and {length_unit_manager(object_cartesian_distance_from_horizon)} from the Outer Horizon."

    #adding space-time velocity
    if bh_type in ["schw", "1"]:
        fd_text = "No Frame Dragging."
    elif bh_type in ["kerr", "2"]:
        fdl_velocity = object_data[7]
        fd_text = f"Frame-dragging coordinate linear speed: {speed_unit_manager(fdl_velocity)}."
    
    sphere_theta = np.linspace(0, np.pi, 20)
    sphere_phi = np.linspace(0, 2 * np.pi, 20)

    sphere_theta, sphere_phi = np.meshgrid(sphere_theta, sphere_phi)

    x = x_center + (outer_ev / 60) * np.sin(sphere_theta) * np.cos(sphere_phi)
    y = y_center + (outer_ev / 60) * np.sin(sphere_theta) * np.sin(sphere_phi)
    z = z_center + (outer_ev / 60) * np.cos(sphere_theta)

    object = ax.plot_surface(x, y, z, color="black", alpha=1, edgecolor=edgecolor, linewidth=0.15)

    x = x_center / distance
    y = y_center / distance
    z = z_center / distance

    #drawing the arrow for the object
    stationary_arrow = ax.plot([x_center, x_center + (x * (outer_ev / 2))], [y_center, y_center + (y * (outer_ev / 2))], [z_center, z_center + (z * (outer_ev / 2))], color="black", linewidth=1)

    #adding the text for the object
    stationary_text = ax.text(x_center + (x * (outer_ev / 2)) , y_center + (y * (outer_ev / 2)), z_center + (z * (outer_ev / 2)), f"{text} \n {fd_text}", 
                              color="white", path_effects=[pe.withStroke(linewidth=1.6, foreground="black")], fontsize=10)

    ax.set_box_aspect([1, 1, 1])
    
    #show the object
    plt.pause(0.5)

    return stationary_arrow, stationary_text, x_center, y_center, z_center, object

def orbiting_object_xp_draw(bh_type, bh_mass, a_star, distance, direction, s_polar_angle, stationary_arrow, stationary_text, x_center, y_center, z_center, object):
    #removing arrow and text
    for arrow in stationary_arrow:
        arrow.remove()

    stationary_text.remove()

    outer_ev = outer_event_horizon_radius(bh_type, bh_mass, a_star, s_polar_angle, False)

    photon_ring = photon_sphere_radius(bh_type, bh_mass, a_star, direction, s_polar_angle, False)

    match bh_type:
        case "schw" | "1":
            color = "#e0115f"
            text = "Co-Rotating Object's Orbit Cartesian Radius at"

            fda_velocity = 0
            fdl_velocity = 0
        case "kerr" | "2":
            match direction:
                case "1":
                    color = "#e0115f"
                    text = "Co-Rotating Object's Orbit Cartesian Radius at"
                case "2":
                    color = "#1160b0"
                    text = "Counter-Rotating Object's Orbit Cartesian Radius at"

            #adding space-time velocity at object's ring
            fda_velocity = frame_dragging_angular_velocity(bh_mass, a_star, distance, gravitational_radius(bh_mass), s_polar_angle, False)
            fdl_velocity = frame_dragging_linear_velocity(bh_mass, a_star, fda_velocity, distance, s_polar_angle)

    #object's radial direction
    r = np.array([x_center, y_center, z_center], dtype=float)
    r = r / np.linalg.norm(r)

    #first orbital direction
    u = np.cross(r, [0, 0, 1])

    if np.linalg.norm(u) < 1e-6:
        u = np.cross(r, [1, 0, 0])

    u = u / np.linalg.norm(u)

    #create the object's ring
    phi = np.linspace(0, 2 * np.pi, 100)
    radius = np.linalg.norm([x_center, y_center, z_center])

    x = radius * (r[0] * np.cos(phi) + u[0] * np.sin(phi))
    y = radius * (r[1] * np.cos(phi) + u[1] * np.sin(phi))
    z = radius * (r[2] * np.cos(phi) + u[2] * np.sin(phi))

    #plot the object's orbit ring
    ax.plot(x, y, z, color=color, linewidth=2)

    x_text = x_center / distance
    y_text = y_center / distance
    z_text = z_center / distance

    #drawing the arrow for object's orbit ring
    ax.plot([x_center, x_center + (x_text * (outer_ev / 2))], [y_center, y_center + (y_text * (outer_ev / 2))], [z_center, z_center + (z_text * (outer_ev / 2))], color=color)

    #adding the text for object's orbit ring
    ax.text(x_center + (x_text * (outer_ev / 2)) , y_center + (y_text * (outer_ev / 2)), z_center + (z_text * (outer_ev / 2)), f"{text} {round(np.degrees(s_polar_angle))}°: {length_unit_manager(cartesian_radius(bh_mass, distance, a_star, s_polar_angle))} from the Black Hole's Center.\n and {length_unit_manager(cartesian_radius(bh_mass, distance - outer_ev, a_star, s_polar_angle))} from the Outer Horizon.\n Frame-dragging coordinate linear speed: {speed_unit_manager(fdl_velocity)}", color=color, path_effects=[pe.withStroke(linewidth=2, foreground="white")], fontsize=10)

    ax.set_box_aspect([1, 1, 1])
    
    #show the object
    plt.pause(0.5)

#velocity calculators
def frame_dragging_angular_velocity(bh_mass, a_star, distance, gravitational_radius, theta, show):
    if a_star != 0:
        #frame dragging orbital angular velocity calculation
        dimensionless_radius = distance / gravitational_radius

        #Σ* and Δ* (theta-dependent kerr metric parameters)
        sigma_star = (dimensionless_radius ** 2) + ((a_star ** 2) * (np.cos(theta) ** 2))
        delta_star = (dimensionless_radius ** 2) - (2 * dimensionless_radius) + (a_star ** 2)

        #A* (theta-dependent denominator)
        A_metric_star = ((dimensionless_radius ** 2) + (a_star ** 2)) ** 2 - ((a_star ** 2) * delta_star * (np.sin(theta) ** 2))

        part1 = (2 * a_star * dimensionless_radius) / A_metric_star
        part2 = c3 / (gravitational_constant * bh_mass)
        l_omega = part1 * part2

        #frame drgging coordinate linear velocity calculation
        fdl_velocity = frame_dragging_linear_velocity(bh_mass, a_star, l_omega, distance, theta) #flag

        if show == True:
            print()
            print(Fore.RED + "Frame Dragging Informations:" + Fore.RESET)
            print("========================================")
            print(Fore.LIGHTGREEN_EX + f"Frame Dragging Orbital Angular Velocity at this Distance: {Fore.LIGHTBLUE_EX} {l_omega:,.6f} rad/s" + Fore.RESET)
            print()
            print(Fore.LIGHTGREEN_EX + f"Frame Dragging Orbital Period: {Fore.LIGHTBLUE_EX} {orbital_period(l_omega, False):,.6f} Second/s" + Fore.RESET)
            print()
            print(Fore.LIGHTGREEN_EX + f"Space_Time Rotates around the Black Hole: {Fore.LIGHTBLUE_EX} {orbital_frequensy(l_omega, False):,.6f} Times per Second!" + Fore.RESET)
            print()
            print(Fore.LIGHTGREEN_EX + f"Frame-Dragging Coordinate Speed at this Radius: {Fore.LIGHTBLUE_EX} {fdl_velocity:,.6f} Meter/s per Second!" + Fore.RESET)
            print()
            print(Fore.LIGHTBLACK_EX + f"{speed_unit_manager(fdl_velocity)}!" + Fore.RESET)
            print()

        return l_omega

def frame_dragging_linear_velocity(bh_mass, a_star, omega, radius, theta):
    a = a_star * gravitational_radius(bh_mass)
    fdl_velocity = omega * np.sqrt((radius ** 2) + (a ** 2)) * np.sin(theta)

    return fdl_velocity
    
def orbital_angular_velocity(bh_mass, a_star, gravitational_radius, distance, direction, show):
    part1 = c3 / (gravitational_constant * bh_mass)
    dimensionless_radius = distance / gravitational_radius
    match a_star:
        case 0:
            omega = sqrt((gravitational_constant * bh_mass) / (distance ** 3))
        case _:
            match direction:
                case "prograde":
                    omega = part1 * (1 / (dimensionless_radius ** (3/2) + a_star))

                case "retrograde":
                    omega = -part1 * (1 / (dimensionless_radius ** (3/2) - a_star))
    if show == True:
        print(Fore.LIGHTGREEN_EX + f"Orbital Angular Velocity: {Fore.LIGHTBLUE_EX} {omega:,.6f} rad/s" + Fore.RESET)
        print()
    return omega

def orbital_period(orbital_angular_velocity, show):
    orbital_period = (2 * pi) / abs(orbital_angular_velocity)
    if show == True:
        print(Fore.LIGHTGREEN_EX + f"Orbital Period: {Fore.LIGHTBLUE_EX} {orbital_period:,.6f} Second/s" + Fore.RESET)
        print()
    return orbital_period

def orbital_frequensy(orbital_angular_velocity, show):
    orbital_frequensy = abs(orbital_angular_velocity) / (2 * pi)
    if show == True:
        print(Fore.LIGHTGREEN_EX + f"Rotation rate around the Black Hole about: {Fore.LIGHTBLUE_EX} {orbital_frequensy:,.6f} Times per Second!" + Fore.RESET)
        print()
    return orbital_frequensy

def coordinate_orbital_velocity(orbital_angular_velocity, distance, show):
    v_coordinate = abs(orbital_angular_velocity) * distance
    if show == True:
        print(Fore.LIGHTGREEN_EX + f"Orbital Velocity (as seen by Distant Observer): {Fore.LIGHTBLUE_EX} {v_coordinate:,.6f} Meter/s per Second!" + Fore.RESET)
        print(Fore.LIGHTBLACK_EX + f"{speed_unit_manager(v_coordinate)}!" + Fore.RESET)
        print()
    return v_coordinate
    
def local_orbital_velocity(a_star, event_horizon_radius, orbital_angular_velocity, frame_dragging_angular_velocity, distance, gravitational_radius, show):
    if distance - event_horizon_radius == 0:
        print(Fore.RED + f"there is NO Real Orbital Velocity at the Edge of the Event Horizon!!!" + Fore.RESET)
        print()
    else:
        dimensionless_radius = distance / gravitational_radius
        if a_star == 0:
            v_local = c * sqrt(1 / (dimensionless_radius - 2))  
        else:
            #Ω*
            omega_rel = orbital_angular_velocity - frame_dragging_angular_velocity
            #Δ*
            delta_star = (dimensionless_radius ** 2) - (2 * dimensionless_radius) + (a_star ** 2)
            #Σ*
            sigma_star = dimensionless_radius ** 2
            #A*
            A_metric_star = ((((dimensionless_radius ** 2) + (a_star ** 2)) ** 2) - ((a_star ** 2) * (delta_star)))
            #α*
            lapse_function = sqrt((delta_star * sigma_star) / A_metric_star)
            #gΦΦ*
            g_phiphi_star = A_metric_star / sigma_star

            v_local = abs(((omega_rel) * (gravitational_radius) * sqrt(g_phiphi_star)) / (lapse_function))

            if show == True:
                print(Fore.LIGHTGREEN_EX + f"Real Orbital Velocity: {Fore.LIGHTBLUE_EX} {v_local:,.6f} Meter/s per Second!" + Fore.RESET)
                print(Fore.LIGHTBLACK_EX + f"{speed_unit_manager(v_local)}!" + Fore.RESET)
                print()

        return v_local

#time
def time_dilation(a_star, near_type, direction, bh_mass, distance, time, theta):

    #schwarschild black hole and stationary
    if a_star == 0 and near_type == "stationary":
        universe_time = time / (sqrt(1 - (2 * gravitational_constant * bh_mass) / (distance * c2)))

    #schwarschild black hole and orbiting
    elif a_star == 0 and near_type == "orbiting":
        if distance >= r_isco(a_star, "", gravitational_radius(bh_mass),"", False):
            print(Fore.LIGHTBLUE_EX + "You're outside of ISCO. Stable Circular Orbit." + Fore.RESET)
            universe_time = time / (sqrt(1 - (3 * gravitational_constant * bh_mass) / (distance * c2)))
        elif photon_sphere_radius(a_star, gravitational_radius(bh_mass), False) <= distance < r_isco(a_star, "", gravitational_radius(bh_mass),"", False):
            print(Fore.LIGHTRED_EX + "You're Between Photon Sphere and ISCO. Circular Orbit Exists but is Unstable!!" + Fore.RESET)
            universe_time = time / (sqrt(1 - (3 * gravitational_constant * bh_mass) / (distance * c2)))
        else:
            return False

    #kerr black hole and stationary
    elif a_star != 0 and near_type == "stationary":
        rho = (distance * c2) / (gravitational_constant * bh_mass)
        if distance >= outer_ergosphere_radius(bh_mass, a_star, theta, False):
            universe_time = time / (sqrt(1 - (2 * rho) / ((rho ** 2) + (a_star ** 2) * (cos(radians(theta))) ** 2) ))
        else:
            return False

    #kerr black hole and orbiting
    elif a_star != 0 and near_type == "orbiting":
        #calculating dimensionless_radius
        #r*
        dimensionless_radius = distance / gravitational_radius(bh_mass)

        #calculating kerr metric parameters(dimensionless)

        #Σ*
        sigma_star = (dimensionless_radius ** 2) + ((a_star ** 2) * (cos(radians(theta)) ** 2))

        #calculating kerr metric components

        #gtt
        time_component = -(1 - ((2 * dimensionless_radius) / sigma_star))

        #gtΦ*
        time_azimuthal_component = -(((2 * a_star) * (dimensionless_radius) * (sin(radians(theta)) ** 2)) / (sigma_star))

        #gΦΦ*
        azimuthal_azimuthal_component = (sin(radians(theta)) ** 2) * ((dimensionless_radius ** 2) + (a_star ** 2) + (((2 * (a_star ** 2)) * (dimensionless_radius) * (sin(radians(theta)) ** 2)) / (sigma_star)))

        #Ω*
        omega_star = orbital_angular_velocity(bh_mass, a_star, gravitational_radius(bh_mass), distance, direction, False) * ((gravitational_constant * bh_mass) / c3)

        ##calculating time dilation 
        try:
            universe_time = time / sqrt(-((time_component) + (2 * time_azimuthal_component * omega_star) + (azimuthal_azimuthal_component * (omega_star ** 2))))
        except ValueError or ZeroDivisionError:
            print(Fore.RED + "This Orbit Is Not Physically Possible. No Timelike Circular Path Exists in This Region!" + Fore.RESET)
            return False

        match direction:
            case "prograde":
                if distance >= r_isco(a_star, direction, gravitational_radius(bh_mass), "", False):
                    print(Fore.LIGHTBLUE_EX + "You're outside of Inner ISCO. Stable Circular Orbit." + Fore.RESET)
                elif photon_sphere_radius(a_star, gravitational_radius(bh_mass), False)[0] <= distance < r_isco(a_star, direction, gravitational_radius(bh_mass),"", False):
                    print(Fore.LIGHTRED_EX + "You're Between Inner Photon Sphere and Inner ISCO. Circular Orbit Exists but is Unstable!!" + Fore.RESET)
                else:
                    print(Fore.RED + "You're inside the Inner Photon Sphere! No Stable Massive Circular Orbit!!! results are Hypothetical." + Fore.RESET)
            case "retrograde":
                if distance >= r_isco(a_star, direction, gravitational_radius(bh_mass), "", False):
                    print(Fore.LIGHTBLUE_EX + "You're outside of Outer ISCO. Stable Circular Orbit." + Fore.RESET)
                elif photon_sphere_radius(a_star, gravitational_radius(bh_mass), False)[0] <= distance < r_isco(a_star, direction, gravitational_radius(bh_mass),"", False):
                    print(Fore.LIGHTRED_EX + "You're Between Outer Photon Sphere and Outer ISCO. Circular Orbit Exists but is Unstable!!" + Fore.RESET)
                else:
                    print(Fore.RED + "You're inside the Outer Photon Sphere! No Stable Massive Circular Orbit!!! results are Hypothetical." + Fore.RESET)

    return universe_time

#calling the engine
if __name__ == "__main__":
    try:
        main()
    except EOFError:
        pass