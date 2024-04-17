
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

def get_compressive_stresses_corners(theta: float, b1:float, t1:float, b2:flot, t2:float, mx:float):
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
        return min_stress
    else:
        return 0

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

if __name__ == "__main__":
    b1 = 100
    t1 = 20
    b2 = 40
    t2 = 30
    I = get_stiffener_inertia(b1, t1, b2, t2)
    print(f"2MOA: {I}")
    print(f"Principle Axes: {np.degrees(get_principal_axes(I[0], I[1], I[2]))}")
    print(f"2MOA after 90 degree rotation: {rotate_I_axes(np.pi/2, I[0], I[1], I[2])}")
    print(f"2MOA after 180 degree rotation: {rotate_I_axes(np.pi, I[0], I[1], I[2])}")


