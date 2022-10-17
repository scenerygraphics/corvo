from PyQt5.QtCore import pyqtSlot
from PyQt5.QtWidgets import QWidget, QMessageBox


class ConfirmPopup(QWidget):
    """
    Generic popup confirmation wrapper around arbitrary function and args.
    """
    def __init__(self, message: str, fn, *args):
        super(ConfirmPopup, self).__init__()

        self.message = message
        self.fn = fn
        self.args = args

        msg = QMessageBox()
        msg.setWindowTitle("Corvo Launcher")
        msg.setText(self.message)
        msg.setIcon(QMessageBox.Question)
        msg.setStandardButtons(QMessageBox.Cancel | QMessageBox.Yes)
        msg.setDefaultButton(QMessageBox.Yes)

        msg.buttonClicked.connect(self.handle_click_event)
        msg.exec()

    def handle_click_event(self, button):
        if button.text() == "&Yes":
            self.fn(*self.args)
        else:
            print("cancelled")
