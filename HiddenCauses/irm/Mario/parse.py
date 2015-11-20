import sys
import re
import random
vevent = re.compile("V\s+([\-\d\w]+)\s+(\w+)")
sevent = re.compile("S\s+([\-\d\w]+)\s+([\-\d\w]+)")
aevent = re.compile("A\s+([\-\d\w]+)")
devent = re.compile("D\s+([\-\d\w]{2,})")
cevent = re.compile("C\s+([\-\d\w]+)\s+[\-\d]*\s*([\-\d\w]+)\s+([\-\d\w]+)")
animate = ["Goomba","Mario","BigMario","FireMario","GreenKoopa","RedKoopa"]
enemies = ["Goomba","GreenKoopa","RedKoopa"]
dirs = ["U","D","L","R"]
opposite = {
	"U":"D",
	"D":"U",
	"L":"R",
	"R":"L"
}
enemyOdds = 1.0/3200.0
bushOdds = 1.0/3200.0
with open(sys.argv[1],'r') as openfile:
	print "ObjectA,ObjectB,A2BDir,EffectType,Source,Target,VelChange"
	causes = []
	effects = []
	for line in openfile:
		if 'NEWFRAME' in line:
			#print causes
			if random.random() < bushOdds:
				an = random.choice(animate)
				d =random.choice(dirs)
				causes.append(["Bush",an,d])
				causes.append([an,"Bush",opposite[d]])
			if random.random() < enemyOdds:
				e1 = random.choice(enemies)
				e2 = random.choice(enemies)
				d =random.choice(dirs)
				causes.append([e1,e2,d])
				causes.append([e2,e1,opposite[d]])
			if not causes:
				pass
				#causes.append(["None","None","None"])
			for cause in causes:
				if not effects:
					print ",".join(cause) + ",None,None,None,None"
				for effect in effects:
					print ",".join(cause) + "," + ",".join(effect)
			causes = []
			effects = []
		else:
			amatch = aevent.search(line)
			dmatch = devent.search(line)
			smatch = sevent.search(line)
			cmatch = cevent.search(line)
			vmatch = vevent.search(line)
			if amatch:
				effects.append(["append",amatch.group(1),"None","None"])
			if vmatch:
				effects.append(["VelChange",vmatch.group(1),"None",vmatch.group(2)])
			if smatch:
				effects.append(["Change",smatch.group(1),smatch.group(2),"None"])
			if dmatch:
                                if 'RUN' not in line:
				        effects.append(["Delete",dmatch.group(1),"None","None"])
			if cmatch:
				o1 = cmatch.group(1)
				o2 = cmatch.group(2)
				if "-" in o1:
					o1 = "B" + o1
				if "-" in o2:
					o2 = "B" + o2
				causes.append([o1,o2,cmatch.group(3)])
				causes.append([o2,o1,opposite[cmatch.group(3)]])



