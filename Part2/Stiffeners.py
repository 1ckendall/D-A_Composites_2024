
# All coordinates defined for an inverted T structure, with Y axis on vertical symmetry plane and X axis on bottom surf.
# 1 is point of T, 2 is base. 2 is 2 * b2 wide!

import numpy as np

def get_rectangle_area(a, b):
    return a*b

def get_total_stiffener_area(b1, t1, b2, t2):
    area_1 = get_rectangle_area(b1, t1)
    area_2 = get_rectangle_area(2 * b2, t2)
    return area_1 + area_2

def get_stiffener_centroid(b1, t1, b2, t2):
    area_1 = get_rectangle_area(b1, t1)
    area_2 = get_rectangle_area(2 * b2, t2)
    return (area_1 * (t2 + 0.5 * b1) + area_2 * 0.5 * t2) / (area_1 + area_2)


def get_Ixx(width, height):
    return (1/12) * width * height**3


def get_Iyy(width, height):
    return (1/12) * height * width**3

def parallel_axis_theorem(I_c, area, d):
    return I_c + area * d**2

def parallel_axis_theorem_Ixy(area, dx, dy):
    return area*dx*dy

def get_stiffener_inertia(b1, t1, b2, t2):
    area_1 = get_rectangle_area(b1, t1)
    area_2 = get_rectangle_area(2*b2, t2)
    centroid_1 = t2 + 0.5 * b1
    centroid_2 = 0.5*t2
    centroid = get_stiffener_centroid(b1, t1, b2, t2)
    y_distance_1 = centroid_1 - centroid
    y_distance_2 = centroid_2 - centroid
    Ixx_1 = get_Ixx(t1, b1)
    Ixx_2 = get_Ixx(2*b2, t2)
    Iyy_1 = get_Iyy(t1, b1)
    Iyy_2 = get_Iyy(2*b2, t2)
    Ixx = parallel_axis_theorem(Ixx_1, area_1, y_distance_1) + parallel_axis_theorem(Ixx_2, area_2, y_distance_2)
    Iyy = Iyy_1 + Iyy_2
    Ixy = parallel_axis_theorem_Ixy(area_1, 0, y_distance_1) + parallel_axis_theorem(area_2, 0, y_distance_2)
    return Ixx, Iyy, Ixy

def rotate_I_axes(theta, Ixx, Iyy, Ixy):
    Ixx_prime = 0.5 * (Ixx + Iyy) + 0.5 * (Ixx - Iyy) * np.cos(2 * theta) - Ixy * np.sin(2 * theta)
    Iyy_prime = 0.5 * (Ixx + Iyy) - 0.5 * (Ixx - Iyy) * np.cos(2 * theta) + Ixy * np.sin(2 * theta)
    Ix_prime_y_prime = 0.5 * (Ixx - Iyy) * np.sin(2 * theta) + Ixy * np.cos(2 * theta)
    return Ixx_prime, Iyy_prime, Ix_prime_y_prime


def get_principal_axes(Ixx, Iyy, Ixy):
    """
    Unnecessary function to check that stuff is working properly. For our shape, the principal axes are along the
    symmetry axis and perpendicular to it through the centroid.
    """
    theta_principle_1 = np.arctan((-2 * Ixy)/(Ixx - Iyy))/2
    theta_principle_2 = theta_principle_1 + np.pi / 2
    return theta_principle_1, theta_principle_2

def get_decomposed_m_x_for_rotated_stringer(theta, mx):
    """
    Decompose an applied moment in GLOBAL x into one that applies to a local stringer rotated theta radians ACW.
    """
    mxprime = mx * np.cos(theta)
    myprime = -mx * np.sin(theta)
    return mxprime, myprime

def get_bending_stress(x, y, Mx, My, Ixx, Iyy):
    sigma = (My * -x) / Ixx + (Mx * y) / Iyy
    return sigma

def get_compressive_stresses_corners(theta: float, b1:float, t1:float, b2:float, t2:float, mx:float):
    Ixx, Iyy, _ = get_stiffener_inertia(b1, t1, b2, t2)
    mxprime, myprime = get_decomposed_m_x_for_rotated_stringer(theta, mx)
    centroid = get_stiffener_centroid(b1, t1, b2, t2)
    point0 = (0, t2 + b1 - centroid)
    point1 = (-b2, 0.5 * t2 - centroid)
    point2 = (b2, 0.5 * t2 - centroid)
    stress0 = get_bending_stress(point0[0], point0[1], mxprime, myprime, Ixx, Iyy)
    stress1 = get_bending_stress(point1[0], point1[1], mxprime, myprime, Ixx, Iyy)
    stress2 = get_bending_stress(point2[0], point2[1], mxprime, myprime, Ixx, Iyy)
    most_compressive_stress = min(stress0, stress1, stress2)
    if most_compressive_stress < 0:
        return most_compressive_stress
    else:
        return 0

def get_filler_area(r, turning_flange_thickness):
    filler_area = 2 * (r + turning_flange_thickness)**2 * (1 - np.pi / 4)
    return filler_area

def NEF_crippling(b, t, sigma_ult_c = 1):
    """B-basis NEF Crippling"""
    if b < 8.443 * t:
        return sigma_ult_c
    else:
        return sigma_ult_c * 11.0 / ((b / t)**1.124)

def OEF_crippling(b, t, sigma_ult_c = 1):
    """B-basic OEF Crippling"""
    if b < 1.98 * t:
        return sigma_ult_c
    else:
        return sigma_ult_c * 1.63 / ((b / t)**0.717)


def column_buckling(K, E, I, L):
    P_cr = -(np.pi**2 * E * I) / ((K * L)**2)
    return P_cr

def column_buckling_elastic_foundation(E, I, L, k):
    column_buckling_ss_P_cr = column_buckling(1, E, I , L)
    m = 1
    elastic_foundation_P_cr = np.inf
    while 1:
        P_cr = column_buckling_ss_P_cr * (m**2 + k * L**4 / (np.pi**4 * E * I * m**2))
        if P_cr < elastic_foundation_P_cr:
            elastic_foundation_P_cr = P_cr
            m+=1
        else:
            break
    return elastic_foundation_P_cr, m


def plate_buckling(b, t, D66):
    Nx_crit = -12 * D66 / b**2
    sigma_crit = Nx_crit / t
    return sigma_crit


def symmetric(l):
    l.extend(l.reverse())


if __name__ == "__main__":
    from Class import Lamina, Laminate

    # Composite (UD tape) material properties
    Ex = 142 * 10 ** 9  # [Pa]
    Ey = 11.2 * 10 ** 9  # [Pa]
    Gxy = 5 * 10 ** 9  # [Pa]
    vxy = 0.3
    Xt = 2200 * 10 ** 6  # [Pa]
    Xc = 1800 * 10 ** 6  # [Pa]
    Yt = 70 * 10 ** 6  # [Pa]
    Yc = 300 * 10 ** 6  # [Pa]
    S = 100 * 10 ** 6  # [Pa]
    rho = 1610  # [kg/m^3]
    t_ply = 0.135 * 10 ** -3  # [m]

    #### Stringer A

    web_layup =     [45, -45, 45, -45, 0, 0, 45, -45, 0, 90, 90, 0, -45, 45, 0, 0, -45, 45, -45, 45]
    flange_layup =  [45, -45, 45, -45, 45, -45, 45, -45, 90, 0,  90, 90, 0, 90, -45, 45, -45, 45, -45, 45, -45, 45]

    # Web
    b1 = 20E-3
    # Flange
    b2 = 20E-3
    # Joint radius
    r = 4E-3
    # Length
    length = 500E-3
    # Effective length factor
    K = 0.5  # 0.5 for fixed ends, 1 for free ends

    # #### Stringer B
    #
    # web_layup = [45,-45,0,90,90,-60,60,-30,30,0,-45,45,45,-45,0,-90,-90,60,-60,30,-30,0,45,-45]
    # flange_layup = [45, -45, 45, -45, -45, 45, 0, 0, 90, 90, 0, 0, 45, -45, 45, -45, 45, -45]
    #
    # # Web
    # b1 = 40E-3
    # # Flange
    # b2 = 40E-3
    # # Joint radius
    # r = 4E-3
    # # # Length
    # length = 500E-3
    # # # Effective length factor
    # K = 0.5  # 0.5 for fixed ends, 1 for free ends

    # ##### Stringer C
    #
    # web_layup = [45, -45, -45, 45]
    # flange_layup = [45, -45, 0, -45, 45]
    #
    # # Web
    # b1 = 30E-3
    # # Flange
    # b2 = 40E-3
    # # Joint radius
    # r = 4E-3

    t1 = len(web_layup) * t_ply
    t2 = len(flange_layup) * t_ply

    print(f"Web thickness {t1}")
    print(f"Flange thickness {t2}")

    I = get_stiffener_inertia(b1, t1, b2, t2)

    web_plylist = [Lamina(i, Ex, Ey, Gxy, vxy, Xt, Xc, Yc, Yt, S, t_ply) for i in web_layup]
    flange_plylist = [Lamina(i, Ex, Ey, Gxy, vxy, Xt, Xc, Yc, Yt, S, t_ply) for i in flange_layup]

    web = Laminate(web_plylist)
    flange = Laminate(flange_plylist)

    web_area = get_rectangle_area(b1, t1)
    flange_area = get_rectangle_area(2*b2, t2)
    filler_area = get_filler_area(r, 2*t_ply)

    web_E = web.Ex
    flange_E = flange.Ex

    web_EA = web_E * web_area
    flange_EA = flange_E * flange_area
    filler_EA = Ex * filler_area
    total_EA = web_EA + flange_EA + filler_EA

    contribution_web = web_EA / total_EA
    contribution_flange = flange_EA / total_EA
    contribution_filler = filler_EA / total_EA

    total_area = web_area + flange_area + 2 * (1 - np.pi / 4) * r**2
    total_E = total_EA / total_area

    print(f"Total EA: {total_EA/1E9} GPa")
    print(f"Total EI: {total_E * I[0]}")
    print(f"Total E: {total_E/1E9} GPa")
    print(f"Total Area: {total_area} m^2")
    print(f"Moments of inertia: {I} m^4")
    # print(f"Contribution web: {contribution_web}")
    # print(f"Contribution flange: {contribution_flange}")
    # print(f"Contribution filler: {contribution_filler}")
    print()

    sigma = 0
    web_fail = ""
    flange_fail = ""
    web_fail_sigma = None
    flange_fail_sigma = None
    while not web_fail or not flange_fail:
        F = sigma * total_area
        F_web = F * contribution_web
        F_flange = F * contribution_flange
        Nx_web = F_web / b1
        Nx_flange = F_flange / (2 * b2)

        iter_web = Laminate(web_plylist, Nx = Nx_web)
        iter_flange = Laminate(flange_plylist, Nx = Nx_flange)

        iter_web.getStressStrain()
        iter_flange.getStressStrain()

        if not web_fail:
            for i, angle in enumerate(web_layup):
                L = Lamina(angle, Ex, Ey, Gxy, vxy, Xt, Xc, Yt, Yc, S, t_ply, sigma_1=iter_web.sigma11[i])
                L.maxStressFibreFail() #fiber failure
                L.maxStressInterFibreFail() #InterFiberfailure
                L.maxStressShearFail() #Shearfailure
                if L.failuremode:
                    web_fail_sigma = sigma
                    web_fail = L.failuremode

        if not flange_fail:
            for i, angle in enumerate(flange_layup):
                L = Lamina(angle, Ex, Ey, Gxy, vxy, Xt, Xc, Yt, Yc, S, t_ply, sigma_1=iter_flange.sigma11[i])
                L.maxStressFibreFail() #fiber failure
                L.maxStressInterFibreFail() #InterFiberfailure
                L.maxStressShearFail() #Shearfailure
                if L.failuremode:
                    flange_fail_sigma = sigma
                    flange_fail = L.failuremode
        sigma -= 1E6

    print(f"Web compression fail at {web_fail_sigma/1E6} MPa, failure mode {web_fail}")
    print(f"Flange compression fail at {flange_fail_sigma/1E6} MPa, failure mode {flange_fail}")
    print()

    P_crip_web = OEF_crippling(b1, t1, web_fail_sigma)
    P_crip_flange = OEF_crippling(b2, t2, flange_fail_sigma)
    P_col_buckle = column_buckling(K, total_E, I[0], length)

    print(f"Web crippling fail at {P_crip_web / 1E6} MPa")
    print(f"Flange crippling fail at {P_crip_flange / 1E6} MPa")
    # print(f"Column buckling fail at {P_col_buckle / 1E6} MPa")
    print()

    P_plate_buckle_web = plate_buckling(b1, t1, web.ABD[5][5])
    P_plate_buckle_flange = plate_buckling(b2, t2, flange.ABD[5][5])

    print(f"Web plate buckling fail at {P_plate_buckle_web / 1E6} MPa")
    print(f"Flange plate buckling fail at {P_plate_buckle_flange / 1E6} MPa")






