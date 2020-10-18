import time


class Timer:
    def __init__(self):
        self.now = time.time()
        self.init = self.now
        self.begin = 0
        self.end = 0

    def start(self):
        self.begin = time.time()

    def stop(self):
        self.time_convert(time.time() - self.begin)

    def lapsed(self):
        print('Total ', end="")
        self.time_convert(time.time() - self.init)

    @staticmethod
    def time_convert(sec):
        mins = sec // 60
        sec = sec % 60
        hours = int(mins // 60)
        mins = int(mins % 60)
        print("Time Lapsed = {0}:{1}:{2}".format(str(hours).zfill(2), str(mins).zfill(2), round(sec, 2)))
