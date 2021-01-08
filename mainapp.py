from PyQt5.QtWidgets import *
import sys
import numpy as np
import matplotlib.pyplot as plt
from itertools import product, combinations
import math
from scipy.spatial.transform import Rotation as R
from matplotlib.backends.backend_qt5agg import (NavigationToolbar2QT)

from Windows import mainwindow
import Resources

from mplwidget import MplWidget

#######################################################################################################


# equation for hkl lambda
def eqn4(a, h, k, l, a11, a12, a13):
    top = np.abs(np.add(np.add(np.multiply(h, a11), np.multiply(k, a12)), np.multiply(l, a13)))
    denom = np.add(np.add(np.power(h, 2), np.power(k, 2)), np.power(l, 2))
    lambda_calc = 2*a*top/denom
    return lambda_calc


def apply_rotation_matrix(matrix, vector):
    # matrix 3x3
    # vector: list[0,1,2] transpose
    new_vector = []
    new_a = np.add(np.add(np.multiply(matrix[0,0], vector[0]), np.multiply(matrix[0,1], vector[1])), np.multiply(matrix[0,2], vector[2]))
    new_b = np.add(np.add(np.multiply(matrix[1,0], vector[0]), np.multiply(matrix[1,1], vector[1])), np.multiply(matrix[1,2], vector[2]))
    new_c = np.add(np.add(np.multiply(matrix[2,0], vector[0]), np.multiply(matrix[2,1], vector[1])), np.multiply(matrix[2,2], vector[2]))
    new_vector.append(new_a)
    new_vector.append(new_b)
    new_vector.append(new_c)
    # new_vector = np.ndarray([new_a, new_b, new_c])        #to match current s, e?
    return new_vector




class MainApp(QMainWindow, mainwindow.Ui_MainWindow):
    def __init__(self, parent=None):
        super(MainApp, self).__init__(parent)
        self.setupUi(self)
        self.setWindowTitle("Bragg Dip Orientation Visualisation")  # you can do this through qt designer

        #self.addToolBar(NavigationToolbar2QT(self.widget_orientation.canvas, self))
        #self.addToolBar(NavigationToolbar2QT(self.widget_graph.canvas, self))

        # please be a loop onegaishimasu
        self.comboBox_material.addItem("Ni")
        self.comboBox_material.addItem("Ag")
        self.comboBox_material.addItem("Al")
        self.comboBox_material.addItem("alpha_Fe")
        self.comboBox_material.addItem("Cd")
        self.comboBox_material.addItem("Cementite")
        self.comboBox_material.addItem("CeO2")
        self.comboBox_material.addItem("Cu")
        self.comboBox_material.addItem("H2O")
        self.comboBox_material.addItem("Mg")
        self.comboBox_material.addItem("Na")
        self.comboBox_material.addItem("NaCl")
        self.comboBox_material.addItem("Pb")
        self.comboBox_material.addItem("Si")
        self.comboBox_material.addItem("Ti")
        self.comboBox_material.addItem("V")
        self.comboBox_material.addItem("Zn")

        ## you can put these things anywhere (THIS IS THE ONLY CHANGE I MADE)
        self.ax_orientation = self.widget_orientation.figure.add_subplot(111, projection='3d')
        self.ax_orientation.plot([1, 2, 3], [4, 5, 6], [7, 8, 9])   # sample 3d plot
        self.ax_graph = self.widget_graph.figure.add_subplot(111)
        self.ax_graph.plot([1, 2, 3], [4, 5, 6]) # sample line plot

        self.pushButton_plot.clicked.connect(self.update_graph)


    def update_graph(self):
        try:
            ##get orientations and material

            # beam orientation
            # r = R.from_quat([0, 0, np.sin(np.pi/4), np.cos(np.pi/4)])
            b_r = R.from_quat([1, 0, 0, 0])
            b_matrix = b_r.as_matrix()

            # sample orientation
            if self.radioButton_quat.isChecked():
                try:
                    quaternion_angle = []
                    quaternion_angle.append(self.doubleSpinBox_q1.value())
                    quaternion_angle.append(self.doubleSpinBox_q2.value())
                    quaternion_angle.append(self.doubleSpinBox_q3.value())
                    quaternion_angle.append(self.doubleSpinBox_q4.value())
                    r = R.from_quat(quaternion_angle)
                    matrix = r.as_matrix()
                except:
                    print("Bad quaternion input")
            elif self.radioButton_euler.isChecked():
                try:
                    euler_angle = []
                    euler_angle.append(self.doubleSpinBox_e1.value())
                    euler_angle.append(self.doubleSpinBox_e2.value())
                    euler_angle.append(self.doubleSpinBox_e3.value())
                    r = R.from_euler("zyx", euler_angle, degrees=True)
                    matrix = r.as_matrix()
                except:
                    print("Bad Euler input")
            else:
                r = R.from_euler("zyx", [0,0,0], degrees=True)
                matrix = r.as_matrix()

            material = self.comboBox_material.currentText() + ".hkl"
            try:
                h, k, l, m, d, f = np.loadtxt("Resources/{}".format(material), dtype=float, skiprows=1, unpack=True)
                lambda_max = np.multiply(2, d)  ##max lambda for all hkl
                lambda_result = eqn4(3.5195, h=h, k=k, l=l, a11=matrix[0, 0], a12=matrix[0, 1],
                                     a13=matrix[0, 2])  ##hkl lambda for orientation
                cos_alpha = np.divide(lambda_result, lambda_max)
                alpha_rad = np.arccos(cos_alpha)
                alpha_deg = np.multiply(alpha_rad, (180 / np.pi))  ##calc angular change between lambda max and calc
                for index, value in enumerate(alpha_deg):
                    if math.isnan(value):
                        print(str(material) + " aligned with {} {} {} direction".format(int(h[index]), int(k[index]), int(l[index])))
                    else:
                        pass
            except Exception as e:
                self.error = "CRITICAL EXCEPTION: {}".format(e)
                print(str(self.error))





            ##2D graph
            x, y = np.loadtxt("spring_dips_clipped_2.txt", dtype=float, skiprows=1, unpack=True)
            self.widget_graph.canvas.axes.clear()
            self.widget_graph.canvas.axes.plot(x, y, "k.", markersize=2)
            for index, value in enumerate(lambda_result):
                self.widget_graph.canvas.axes.axvline(x=value, color="r")
            for index, value in enumerate(lambda_max):
                self.widget_graph.canvas.axes.axvline(x=value, color="g")
            self.widget_graph.canvas.draw()

            ##3D graph
            # h, k, l, m, d, f = np.loadtxt("{}.hkl".format(self.comboBox_material.activated[str]), dtype=float, skiprows=1, unpack=True)
            #
            # self.widget_orientation.canvas.axes.clear()
            # self.widget_orientation.canvas.axes.plot([-4,4],[0,0],zs=[0,0], color="black", markersize=5)
            # self.widget_orientation.canvas.draw()


        except MemoryError as e:
            self.error = "Unable to calculate strain map: {}".format(e)
        except KeyboardInterrupt:
            raise
        except Exception as e:
            self.error = "CRITICAL EXCEPTION: {}".format(e)
        # finally:
        #     self.finished.emit("Failed" if self.error is not None else ("Stopped" if self.exiting else "Completed"))




def main():
    app = QApplication(sys.argv)
    form = MainApp()
    form.show()
    app.exec_()


if __name__ == '__main__':
    main()
