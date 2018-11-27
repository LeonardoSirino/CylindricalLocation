import kivy
import time
import threading
import numpy as np
from kivy.app import App
from kivy.clock import Clock
from kivy.uix.popup import Popup
from kivy.uix.progressbar import ProgressBar
from kivy.uix.widget import Widget
from kivy.properties import ObjectProperty
kivy.require("1.9.1")

progress = 0
exitFlag = False
t = 0
dt = .005


class MyPopupProgressBar(Widget):

    # Kivy properties classes are used when you create an EventDispatcher.
    progress_bar = ObjectProperty()

    def __init__(self, **kwa):
        # super combines and initializes two widgets Popup and ProgressBar
        super(MyPopupProgressBar, self).__init__(**kwa)
        self.progress_bar = ProgressBar()  # instance of ProgressBar created.
        # progress bar assigned to popup
        self.popup = Popup(title='Barra de progresso',
                           content=self.progress_bar)
        self.popup.bind(on_open=self.puopen)  # Binds super widget to on_open.
        # Uses clock to call progress_bar_start() (callback) one time only
        Clock.schedule_once(self.progress_bar_start)

    # Provides initial value of of progress bar and lanches popup
    def progress_bar_start(self, instance):
        self.progress_bar.value = 1  # Initial value of progress_bar
        self.popup.open()  # starts puopen()

    def next(self, dt):  # Updates Project Bar
        global exitFlag
        global progress

        if exitFlag:  # Checks to see if progress_bar.value has met 100
            return False  # Returning False schedule is canceled and won't repeat

        self.progress_bar.value = progress  # Updates progress_bar's progress

    def puopen(self, instance):  # Called from bind.
        # Creates Clock event scheduling next() every 5-1000th of a second.
        Clock.schedule_interval(self.next, .0005)


class MyApp(App):
    def build(self):
        return MyPopupProgressBar()


def Data():
    global t
    global dt
    global progress
    global exitFlag

    while not exitFlag:
        t += dt
        progress = -50 * np.cos(2 * np.pi * t) + 50
        time.sleep(dt)
        if t > 5:
            exitFlag = True
            return False


def Gui():
    MyApp().run()


t1 = threading.Thread(target=Data)
t2 = threading.Thread(target=Gui)

t1.start()
t2.start()