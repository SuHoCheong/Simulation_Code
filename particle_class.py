class particle:

    def __init__(self,x0,y0,theta):
        self.x0 = x0
        self.y0 = y0
        self.theta = theta
        self.newx0 = x0
        self.newy0 = y0
        self.newtheta = theta

    def setxcoords(self,x0):
        self.x0 = x0 

    def setycoords(self, y0):
        self.y0 = y0

    def getxcoords(self):
        return self.x0

    def getycoords(self):
        return self.y0

    def settheta(self, theta):
        self.theta = theta

    def gettheta(self):
        return self.theta

    def getNEWy(self):
        return self.newy0

    def getNEWx(self):
        return self.newx0

    def getNEWtheta(self):
        return self.newtheta

    def setNEWy(self, y0):
        self.newy0 = y0

    def setNEWx(self, x0):
        self.newx0 = x0

    def setNEWtheta(self, theta):
        self.newtheta = theta
