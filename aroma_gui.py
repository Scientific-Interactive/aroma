#! /usr/bin/env python
# Author : Anuja
# Last Updated : 02.11.2014
# Python 3 upgrade: V. Ganesh
# Updated: 19.10.2020
#

import aroma
import aroma_su

import tkinter as tk
import tkinter.filedialog as tkFileDialog
import tkinter.scrolledtext as ScrolledText
import PIL
from PIL import Image, ImageTk
import threading
import sys
import time

class RedirectText(object):
    def __init__(self, text_ctrl):
        self.output = text_ctrl

    def write(self, string):
        self.output.insert(tk.END, string)


class Application(tk.Frame):
    def __init__(self, master):
        self.colour = "white"
        tk.Frame.__init__(self, master, bg=self.colour)
        self.armFilePath = ""
        self.statusText = tk.StringVar()
        self.statusText.set("Idle") 
        self.grid()
        self.arrangeStaticLabels()
        self.getTheArmFile()
        self.RunAroma()
        self.RunAromaSu()
        self.outputText()
        self.status()
        self.clear()
        self.quitbutton()


    def arrangeStaticLabels(self):

#       Aroma Logo
        self.aroma_img = ImageTk.PhotoImage(   (Image.open("aroma_logo.png")).resize((450,125) ,Image.ANTIALIAS) )
        self.aroma_logo = tk.Label(self, image=self.aroma_img, bg=self.colour)
        self.aroma_logo.image = self.aroma_img
        self.aroma_logo.grid(column=0, row=0, columnspan=3, sticky='NEWS', pady=20 )

        self.inputlabel = tk.Label(self, text="Aroma Input: ", font=("papyrus",12), bg=self.colour)
        self.inputlabel.grid(column=0,row=1,padx=10, pady=10)

        self.statuslabel = tk.Label(self, text="Status: ", font="papyrus", bg=self.colour)
        self.statuslabel.grid(column=0,row=5, padx=10, pady=10, sticky="ew")

#       Technion Logo
        self.tech_img = ImageTk.PhotoImage(   (Image.open("tech_logo.png")).resize((130,100),Image.ANTIALIAS) )
        self.tech_logo = tk.Label(self, image=self.tech_img, bg=self.colour)
        self.tech_logo.image = self.tech_img
        self.tech_logo.grid(column=0, row=6, sticky="news", padx=10, pady=20)

#       Chemistry Logo
        self.chem_img = ImageTk.PhotoImage(   (Image.open("chem_logo.png")).resize((250,100),Image.ANTIALIAS) )
        self.chem_logo = tk.Label(self, image=self.chem_img, bg=self.colour)
        self.chem_logo.image = self.chem_img
        self.chem_logo.grid(column=2, row=6, sticky="news", padx=10, pady=20)

#       Credits
        self.credit = tk.Label(self, text="Citation:\n1. Aroma, Anuja P. Rahalkar and Amnon Stanger\n2. A. Stanger, J. Org. Chem. 2010, 71, 883-893\n3. A. Stanger, J. Org. Chem. 2010, 75, 2281-2288\n4.R. Gershoni-Poranne, A. Stanger, Chem. Eur. J. 2014, 20, 5673-5688", bg="white", justify="center", font=("papyrus", 10))
        self.credit.grid(column=1, row=6, padx=10, pady=10, sticky="ew")

    def getTheArmFile(self):
        self.inpfl = tk.Entry(self)
        self.inpfl.grid(column=1,row=1, sticky="ew")
        self.inpflon = tk.Button(self, text='Browse', command=self.selectFile)
        self.inpflon.grid(column=2,row=1)

    def selectFile(self):
        self.armfl = tkFileDialog.askopenfilename()
        self.inpfl.insert(0,self.armfl)

    def RunAroma(self):
        self.RunSButton = tk.Button(self, text='Run Single', command=self.threadRunS)
        self.RunSButton.grid(column=0,row=2)

    def threadRunS(self):
        aromathread = threading.Thread(target=self.RunS, args=())
        aromathread.start()

    def RunS(self):
        self.fl = tk.Entry.get(self.inpfl)
        self.RunSButton['state']='disabled'
        self.RunMButton['state']='disabled'
        self.statusText.set("Currently Running: " + self.fl)
        aroma.aroma(self.fl[0:len(self.fl)-4])
        self.RunSButton['state']='normal'
        self.RunMButton['state']='normal'
        self.statusText.set("Job Over: " + self.fl)

    def RunAromaSu(self):
        self.RunMButton = tk.Button(self, text='Run Multiple', command=self.threadRunM)
        self.RunMButton.grid(column=1,row=2)

    def threadRunM(self):
        suaromathread = threading.Thread(target=self.RunS, args=())
        suaromathread.start()

    def RunM(self):
        self.fl = tk.Entry.get(self.inpfl)
        self.RunSButton['state']='disabled'
        self.RunMButton['state']='disabled'
        self.statusText.set("Currently Running: " + self.fl)
        aroma_su.main(self.fl[0:len(self.fl)-4])
        self.RunSButton['state']='normal'
        self.RunMButton['state']='normal'
        self.statusText.set("Job Over: " + self.fl)

    def outputText(self):
        self.optext = ScrolledText.ScrolledText(self, height=12, width=100 )
        self.optext.grid(column=0, row=4, columnspan=3, rowspan=1, padx=20, pady=5)
        sys.stdout = RedirectText(self.optext)
        sys.stderr = RedirectText(self.optext)
        
    def status(self):
        self.sts = tk.Label(self, textvariable=self.statusText, bg=self.colour)
        self.sts.grid(column=1,row=5, padx=10, pady=10, sticky="ew")
        self.flash()

    def flash(self):
        self.sts.configure(foreground="black")
        self.after(250, self.flash)
        self.sts.configure(foreground="red4")

    def updateStatus(self):
        self.statusText.set("Idle")
        self.sts.config(textvariable=self.statusText)

    def clear(self):
        self.clearButton = tk.Button(self, text='Clear Log', command=self.cleartext)
        self.clearButton.grid(column=2,row=5)

    def cleartext(self):
        self.optext.delete('1.0', tk.END)

    def quitbutton(self):
        self.quitButton = tk.Button(self, text='Quit', command=self.kill)
        self.quitButton.grid(column=2,row=2, padx=10, pady=10)

    def kill(self):
       raise SystemExit

top = tk.Tk()
top.columnconfigure(0, weight=1)
top.rowconfigure(0, weight=1) 
app = Application(top)
app.master.title('AROMA')    
app.mainloop()                  
