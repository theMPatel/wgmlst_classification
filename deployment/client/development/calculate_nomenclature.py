#//$autorun;event=CreateMainWin

#prototype script for calculating strain nomenclature within the BioNumerics database
#use at your own risk
#written by Hannes Pouseele, hannes_pouseele@applied-maths.com
#this code and its intellectual property is owned by Applied Maths NV
#please do not distribute without explicit consent of the owner


import bns
import os
import gzip
import json
import StringIO 
import base64
import numpy as np
from collections import defaultdict
from wgMLST_Client.wgMLSTSchema import Schema
import xml.etree.cElementTree as ET
	
Dlg = bns.Windows.XmlDlgBuilder

class CalcStrainNamesDlg(Dlg.Dialogs):
	def __init__(self, charExperType, nameFieldId, historyFieldId, charViewId, thresholds, splitNames, minPresenceThreshold):
		
		if len(bns.Database.Db.Fields)<=1:
			raise RuntimeError("Please create database fields to hold the strain names and strain names history first!")
		
		Dlg.Dialogs.__init__(self, 'CalcStrainNamesDlg')
		
		cst = bns.Characters.CharSetType(charExperType.Name)		
		
		# member variables to hold the results
		self.selectedField = nameFieldId if nameFieldId in [f.ID for f in bns.Database.Db.Fields] else bns.Database.Db.Fields[0].ID 
		self.selectedHistoryField = historyFieldId if historyFieldId in [f.ID for f in bns.Database.Db.Fields] else bns.Database.Db.Fields[1].ID
		self.selectedCharView = charViewId if charViewId in cst.ViewGetList() else '__All__'
		self.selectedThresholds = thresholds
		self.selectedSplitNames = splitNames
		self.minPresenceThreshold = minPresenceThreshold
		
		# the values to show in the list control
		fieldItems = [Dlg.ListItem(f.DispName, f.ID) for f in bns.Database.Db.Fields]		
		charViewItems = [Dlg.ListItem(cst.ViewGet(f).GetInfo().get('DisplayName', f), f) for f in cst.ViewGetList()]		
		thresholdItems = [Dlg.ListItem(str(float(f)) +'%') for f in thresholds]
		
		# the controls
		self.fieldCtrl = Dlg.Drop('fieldCtrl', fieldItems, 10, 30, canEdit=False)
		self.historyFieldCtrl = Dlg.Drop('historyFieldCtrl', fieldItems, 10, 30, canEdit=False)
		self.charViewCtrl = Dlg.Drop('charViewCtrl', charViewItems, 10, 30, canEdit=False)
		self.thresholdsCtrl = Dlg.SimpleList('thresholdsCtrl', thresholdItems, 5, 10, True)
		self.newThresholdCtrl = Dlg.Input('newThresholdCtrl', 5)
		self.splitNamesCtrl = Dlg.Check('splitNamesCtrl', "Split names over several fields")
		self.minPresenceThresholdCtrl = Dlg.Input('minPresenceThresholdCtrl', 5, default = 95)
		
		self.addButton = Dlg.Button("action", "Add", contextCall = Dlg.ContextCall(self.OnAddThreshold))
		self.removeButton = Dlg.Button("action", "Remove selected thresholds", contextCall = Dlg.ContextCall(self.OnRemoveThreshold))
		
		
		# now define the dialog layout
		grid = [["Subschema:", self.charViewCtrl],
					 ["Distance thresholds (in %):", Dlg.Cell([[self.thresholdsCtrl, Dlg.Cell([[Dlg.Cell([["Threshold:", self.newThresholdCtrl, "%", self.addButton]])], [self.removeButton]])]])],
				["Name field:", self.fieldCtrl], 
				["History field:", self.historyFieldCtrl], 
				['', self.splitNamesCtrl],
				['Minimum presence threshold:', Dlg.Cell([[self.minPresenceThresholdCtrl, '%']])]
		]
		
		simple = Dlg.SimpleDialog(grid, onStart=self.onStart, onOk=self.onOk)
		self.AddSimpleDialog("Calculate strain nomenclature", simple)
		
	def _ThresholdsFromCtrl(self, selectedOnly = False):
		if selectedOnly:
			return [float(l[:-1]) for l in self.thresholdsCtrl.Multi]
		else:
			return [float(l.ID[:-1]) for l in self.thresholdsCtrl.Choices().ItemsList]

	def _ThresholdsToCtrl(self, thresholds):
		self.thresholdsCtrl.SetItems([Dlg.ListItem(str(f) +'%') for f in thresholds])	
		
	def OnAddThreshold(self, args):
		try:
			newThreshold = float(self.newThresholdCtrl.GetValue()) 			
		except:
			raise RuntimeError("Please provide a number as a threshold.")
			
		if newThreshold<0 or newThreshold>100:
			raise RuntimeError("Please provide a number between 0 and 100 as a threshold.")
		
		oldThresholds = set(self._ThresholdsFromCtrl())
		oldThresholds.add(newThreshold)
		newThresholds = sorted([f for f in oldThresholds])
		self._ThresholdsToCtrl(newThresholds)
		
	
	def OnRemoveThreshold(self, args):
		oldThresholds = self._ThresholdsFromCtrl()
		toRemoveThresholds = set(self._ThresholdsFromCtrl(True))
		newThresholds = sorted([f for f in oldThresholds if f not in toRemoveThresholds])
		self._ThresholdsToCtrl(newThresholds)
	
	def onStart(self, args):
		"""Set the selected items"""
		self.fieldCtrl.SetValue(self.selectedField or "")
		self.historyFieldCtrl.SetValue(self.selectedHistoryField or "")
		self.charViewCtrl.SetValue(self.selectedCharView or "")
		self._ThresholdsToCtrl(self.selectedThresholds)
		self.splitNamesCtrl.Checked = self.selectedSplitNames
		self.minPresenceThresholdCtrl.SetValue(str(self.minPresenceThreshold*100.0))
	
	def onOk(self, args):
		"""Get the selected items"""
		self.selectedField = self.fieldCtrl.GetValue()
		self.selectedHistoryField = self.historyFieldCtrl.GetValue()
		self.selectedCharView = self.charViewCtrl.GetValue()
		self.selectedThresholds = self._ThresholdsFromCtrl()
		self.selectedSplitNames = self.splitNamesCtrl.Checked
		self.minPresenceThreshold = float(self.minPresenceThresholdCtrl.GetValue()) / 100.0


#========== DISTANCE FUNCTION ==================================================
def GetDistance(p1, p2):
	common = np.multiply(p1>0, p2>0)
	nCommon = np.sum(common)
	nSame = np.sum(p1[common]==p2[common])
	if nCommon:
		return 100.0 * (float(nCommon) - float(nSame)) / float(nCommon) #nCommon-nSame
	else:
		return 0.0




#========== NAMES OBJECT =======================================================
class NameEvent(object):
	def __init__(self, oldname, newname, key, info):
		self._oldName = oldname
		self._newName = newname
		self._key = key
		self._info = info


class NameHistory(object):
	def __init__(self):
		self._events = []
		pass

	def AddEvent(self, oldName, newName, key, info):
		self._events.append(NameEvent(oldName, newName, key, info))



class Names(object):
	def __init__(self):
		self._names = defaultdict(list)
		self._history = defaultdict(NameHistory)
		self._allnames = defaultdict(list)
		pass
		
	def _NameToStr(self, parts):
		return '.'.join(str(p) for p in parts)

	def _NameFromStr(self, name):
		return [int(p) for p in name.split('.')]
	
	def Save(self, oFile):
		json.dump({
			'names': {k: self._NameToStr(v) for k, v in self._names.iteritems()}, 
			'allnames': self._allnames.keys(), 
			}, oFile)
		
	def Load(self, iFile):
		data = json.load(iFile)
		self._allnames = defaultdict(list)
		for k in data['allnames']:
			self._allnames[k] = self._NameFromStr(k)
		self._names = defaultdict(list)
		for k, v in data['names'].iteritems():
			self._names[k] = self._NameFromStr(v)

	@property
	def History(self):
		return self._history

	def HasName(self, key):
		return key in self._names
		
	def HasResolvedName(self, key):
		return self.HasName(key) and -1 not in self.GetName(key)

	def DropName(self, key):
		if not self.HasName(key): return
		assert(not self.HasResolvedName(key))
		del self._names[key]

	def GetName(self, key):
		return self._names.get(key, [])

	def GetStrName(self, key):
		return '.'.join(str(n) for n in self.GetName(key))

	def GetPart(self, key, i):
		if i<0: return 1
		else:
			if len(self.GetName(key))<=i:
				kkkkk=0
			return self.GetName(key)[i]

	def SetName(self, key, names):
		self._names[key] = names

	def AddToName(self, key, part):
		self._names[key].append(part)
		
	def MakeUndefined(self, key, length):
		self._names[key] = self._names[key] + [-1] * (length-len(self._names[key]))

	def GetNextIdx(self, parts):
		idx = [p[len(parts)] for k, p in self._allnames.iteritems() if p[:len(parts)]==parts and len(p)>len(parts)]
		if len(idx)==0:
			return 1
		return max(idx)+1 if len(idx)>0 else 1

	def IsSame(self, other):
		nameViolations = []
		nameMap = {}

		for key, name in self._names.iteritems():
			#the real other name of the sample
			otherName = other.GetName(key)

			#the predicted other name of the sample (if any)
			strName = '.'.join(str(n) for n in name)
			strOtherName = '.'.join(str(n) for n in otherName)
			mappedName = nameMap.get(strName, None)
			if mappedName is None: nameMap[strName] = otherName
			else:
				if mappedName!=otherName:
					#if predicted name and other name do not match, we have a problem
					nameViolations.append([strName, strOtherName])

		return len(nameViolations) == 0, nameViolations

	def Merge(self, clusterNames, key, info):

		def NameToStr(parts):
			return '.'.join(str(p) for p in parts)

		#find the new name
		"""
		clusterSizes = {}
		for clusterName in clusterNames:
			n = NameToStr(clusterName)
			for key, name in self._names.iteritems():
				if name[:len(clusterName)] == clusterName:
					clusterSizes[n] = clusterSizes.get(n, 0)+1
		"""

		newName = clusterNames[0]
		#maxSize = clusterSizes[NameToStr(clusterNames[0])]
		"""
		for clusterName in clusterNames:
			n = NameToStr(clusterName)
			s = clusterSizes[n]
			if s > maxSize:
				maxSize = s
				newName = clusterName
		"""



		nextIdx  = {}
		cnt = 0
		changes = {}
		for clusterName in clusterNames:
			for key, name in self._names.iteritems():

				if len(name)==len(clusterName):
					if name==clusterName:
						myOldName = NameToStr(self._names[key])
						self._names[key][-1] = newName[-1]
						myNewName = NameToStr(self._names[key])
						if myOldName != NameToStr(newName):
							changes[myOldName] = myNewName

				elif name[:len(clusterName)] == clusterName:
					if clusterName==newName:
						tail = name[len(clusterName)]
						nextIdx[NameToStr(name[:len(newName)+1])]=tail

					myOldName = NameToStr(self._names[key][:len(clusterName)])

					idx = nextIdx.get(NameToStr(name[:len(newName)+1]), cnt+1)
					if idx==tail and clusterName!=newName: idx+=1

					cnt = max(cnt, idx)
					nextIdx[NameToStr(name[:len(newName)+1])] = idx

					self._names[key][len(newName)-1] = newName[-1]
					self._names[key][len(newName)] = idx

					myNewName = NameToStr(self._names[key][:len(clusterName)])
					if myOldName != NameToStr(newName):
						changes[myOldName] = myNewName

		for oN, nN in changes.iteritems():
			if oN!=nN:
				self._history[oN].AddEvent(oN, nN, key, info)

		return newName[-1]
		
	def Finalize(self, key):
		self._allnames[self.GetStrName(key)] = self.GetName(key) 
	


#========== NAMING FUNCTION ====================================================
def CalcName(namedEntries, existingNames, unNamedEntry, distances, thresholds, qcStatus):
	"""
	:param namedEntries: keys you want to use to identify the samples, can be anything but should be
					in the same order as the list of distances
	:param existingNames: existing names and therefore clusters
	:param unNamedEntry: key for the entry to be named
	:param distances: list of distances between the unNamedEntry and the namedEntries
	:param thresholds: list of thresholds
	:return: a Names object holding the hierarchical names of all the samples
	"""

	#create the object that will hold the names
	names = existingNames

	#create the initial cluster: all entries lumped together
	nEntries = len(namedEntries)+1

	#make sure the thesholds are sorted as biggest-first
	myThresholds = thresholds
	myThresholds.sort(key=lambda x: -x)

	startCluster = range(len(namedEntries)+1)
	cluster = startCluster

	for partNr, threshold in enumerate(myThresholds):

		#create the clusters at this level
		subClusters = defaultdict(list)
		for j in cluster:
			if j==len(namedEntries): continue
			assert(names.GetPart(namedEntries[j], partNr)>0)
			subClusters[names.GetPart(namedEntries[j], partNr)].append(j)

		#find the clusters that are close enough
		closestClusters = set()
		for j, cluster in subClusters.iteritems():
			for k in cluster:
				if distances[k]<=threshold:
					closestClusters.add(j)
		
		#assigned to just one cluster, nothing special happening here
		if len(closestClusters)==1:
			theOne = closestClusters.pop()
			cluster = subClusters[theOne]
			names.AddToName(unNamedEntry, theOne)
		#bad quality, jump out 
		elif not qcStatus:
			names.MakeUndefined(unNamedEntry, len(myThresholds))
			break
		#create new cluster
		elif len(closestClusters)==0:
			cluster = []
			idx = names.GetNextIdx(names.GetName(unNamedEntry))
			names.AddToName(unNamedEntry, idx)
		#assigned to multiple clusters, need to merge these clusters	
		else:
			#assigned to multiple clusters, need to merge these clusters
			cluster = []
			clusterSizes = []
			for c in closestClusters:
				cluster.extend(subClusters[c])
				clusterSizes.append([c, len(subClusters[c])])

			clusterSizes.sort(key=lambda x: -x[1])
			prefix = names.GetName(unNamedEntry)
			newName = names.Merge([list(prefix) + [c[0]] for c in clusterSizes], unNamedEntry, None)
			names.AddToName(unNamedEntry, newName)
	
	if qcStatus:
		names.Finalize(unNamedEntry)
	
	return names

class Profiles(object):
	def __init__(self):
		self._profiles = None
		self._hashToIndex = {}
		self._entryToHash = {}
		
	def Load(self, iFile):
		data = json.load(iFile)
		self._entryToHash = data['entryToHash']
		self._hashToIndex = data['hashToIndex']
		iFl =  StringIO.StringIO(base64.b64decode(data['profiles']))
		
		self._profiles = np.load(iFl)
		self._profiles = self._profiles['arr_0']
		
	def Save(self, oFile):
		oFl = StringIO.StringIO()
		np.savez_compressed(oFl, self._profiles)
		profilesStr = oFl.getvalue()
		oFl.close()
		
		json.dump({
			'entryToHash': self._entryToHash,
			'hashToIndex': self._hashToIndex,
			'profiles': base64.b64encode(profilesStr)
			}, oFile)
			
		i=0
	
	def Add(self, h, profile):
		if self._profiles is None:
			self._profiles = profile
		else:
			self._profiles = np.concatenate((self._profiles,  profile), axis = 0)
		self._hashToIndex[h] = len(self._profiles)-1
			
	def Get(self, h):
		return self._profiles[self._hashToIndex[h]]
	
	@property
	def EntryToHash(self):
		return self._entryToHash
	
	
#========== MAIN ==============================================================

class Calculator(object):
	def __init__(self, experType, nameFieldId, nameHistoryFieldId, charViewId, thresholds, splitNames, minPresenceThreshold, allNamesFlName, allProfilesFlName):
		self.experType = experType
		self.nameFieldId = nameFieldId
		self.nameHistoryFieldId = nameHistoryFieldId
		self.charViewId = charViewId
		self.thresholds = thresholds
		self.splitNames = splitNames
		self.minPresenceThreshold = minPresenceThreshold
		
		self.allNamesFlName = allNamesFlName
		self.allProfilesFlName = allProfilesFlName
		
		self.names = None
		
	def DoRefresh(self, args):
		pass
		
	def DoCalc(self, args):
		comm = args.get('communication', bns.Windows.CalcCommunication())
	
		#fetch the wgMLST experiment type
		experType = self.experType
		
		comm.SetMessage("Fetching character subset ...")	
		#fetch the character numbers we need
		charNrs = []
		cst = bns.Characters.CharSetType(experType.Name)
		if self.charViewId is None or len(self.charViewId)==0 or self.charViewId not in cst.ViewGetList():
			charNrs = range(cst.GetCount())
		else:
			charView = cst.ViewGet(self.charViewId)
			chars = set(charView.GetCharacters())
			charNrs = [i for i in range(cst.GetCount()) if cst.GetChar(i) in chars]
		
		#fetch the existing names
		comm.SetMessage("Fetching existing names ...")		
		existingNames = Names()
		if os.path.exists(self.allNamesFlName):
			with gzip.open(self.allNamesFlName, 'rb') as iFile:
				existingNames.Load(iFile)
			
		#fetch the profiles 
		profiles = Profiles()
		if os.path.exists(self.allProfilesFlName):
			with gzip.open(self.allProfilesFlName, 'rb') as iFile:
				profiles.Load(iFile)
		
		entryToHash = profiles.EntryToHash
		namedHashes = [h for k, h in entryToHash.iteritems() if existingNames.HasResolvedName(h)]
		namedHashesSet = set(namedHashes)
		namedEntries = [k for k,h in entryToHash.iteritems() if existingNames.HasResolvedName(h)]
		
		noExperEntries = set()
		badQualityEntries = set()
		
		#fetch the existing characters for the selected entries
		comm.SetMessage("Fetching profiles ...")	
		
		nSelected = len(bns.Database.Db.Selection)
		for i, e in enumerate(bns.Database.Db.Selection):
			if i%100==0: comm.SetProgress(i, nSelected)
			exper = bns.Database.Experiment(e, experType.Name)
			if not exper.IsPresent():
				noExperEntries.add(e.Key)
			else:
				vals = []
				presences = []
				exper.LoadOrCreate().GetData(vals, presences)
				nPresent = sum(presences[i] for i in charNrs)
				if float(nPresent) / float(len(charNrs))<self.minPresenceThreshold:
					badQualityEntries.add(e.Key)
				
				profile = np.asarray([[vals[i] for i in charNrs]])
				profile.flags.writeable = False
				h=str(hash(profile.data))
				
				profiles.Add(h, profile)
				entryToHash[e.Key] = h
				if existingNames.HasName(h):
					namedEntries.append(e.Key)	
					
		#for each entry, calculate the name	
		comm.SetMessage("Calculating names ...")	
		lenSel = len(bns.Database.Db.Selection)
		for i, entry in enumerate(bns.Database.Db.Selection):
			comm.SetProgress(i, lenSel)
	
			qcStatus = entry.Key not in badQualityEntries
				
			if qcStatus and entryToHash[entry.Key] in namedHashesSet:
				continue
			
			if entry.Key in noExperEntries:
				continue
				
			if existingNames.HasName(entryToHash[entry.Key]):
				existingNames.DropName(entryToHash[entry.Key])
		
			myChars = profiles.Get(entryToHash[entry.Key])
			
			#calculate the distance between the unnamed sample and all the named samples
			dists = [GetDistance(myChars, profiles.Get(h)) for h in namedHashes]
			
			#calculate the name of the entry
			existingNames = CalcName(namedHashes, existingNames, entryToHash[entry.Key], dists, self.thresholds, qcStatus)
			
			#keep track of the data
			if qcStatus:
				namedHashesSet.add(entryToHash[entry.Key])
				namedHashes.append(entryToHash[entry.Key])
			
		#store the allnames object
		with gzip.open(self.allNamesFlName, 'wb') as oFile:
			existingNames.Save(oFile)	
			
		with gzip.open(self.allProfilesFlName, 'wb') as iFile:
			profiles.Save(iFile)
		
		
		#save the name of the entry
		comm.SetMessage("Saving names ...")		
		if self.splitNames:
			existingFields = { f.ID: f.DispName for f in bns.Database.Db.Fields }		
			fieldIds = { self.nameFieldId + '_' + str(threshold).replace('.', '_') : existingFields[self.nameFieldId] + ' ({}%)'.format(threshold) for threshold in self.thresholds }
			for fieldId, fieldName in fieldIds.iteritems():
				if fieldId not in existingFields:
					bns.Database.Db.Fields.Add(fieldId, fieldname =  fieldName)
				
		fieldIds = [ self.nameFieldId + '_' + str(threshold).replace('.', '_') for  threshold in self.thresholds ]
		for key, h in entryToHash.iteritems():
			oldName = bns.Database.EntryField(key, self.nameFieldId).Content 
			newName = existingNames.GetStrName(entryToHash[key])
			if oldName==newName: continue
			
			#from now on, oldName!=newName
			if oldName:
				nameHistory = bns.Database.EntryField(key, self.nameHistoryFieldId).Content
				bns.Database.EntryField(key, self.nameHistoryFieldId).Content = oldName + ', ' + nameHistory 
			
			bns.Database.EntryField(key, self.nameFieldId).Content = newName
			
			if self.splitNames:
				name = existingNames.GetName(entryToHash[key])
				for i, fieldId in enumerate(fieldIds):
					bns.Database.EntryField(key, fieldId).Content = str(name[i])
			
		bns.Database.Db.Fields.Save()
	
	
def Main(args):
	
	winId = 1
	
	#pick up previously used settings
	tmpNameFieldId = ''
	tmpNameHistoryFieldId = ''
	tmpCharViewId = ''
	tmpThresholds = [0.5, 1, 5, 10]
	tmpSplitNames = False
	tmpMinPresenceThreshold = 0.95

	settingsKey = "hierarchical_strain_nomenclature"
	settings = bns.Database.Db.Info.LoadSetting(settingsKey, True)
	if settings:
		root = ET.fromstring(settings)
		n = root.find('nameFieldId')
		if n is not None: tmpNameFieldId = n.text
		
		n = root.find('nameHistoryFieldId')
		if n is not None: tmpNameHistoryFieldId = n.text
		
		n = root.find('charViewId')
		if n is not None: tmpCharViewId = n.text
		
		n = root.find('thresholds')
		if n is not None: tmpThresholds = [float(f) for f in n.text.split('|')]

		n = root.find('splitNames')
		if n is not None: tmpSplitNames = n.text == '1'
		
		n = root.find('minPresenceThreshold')
		if n is not None: tmpMinPresenceThreshold = float(n.text)
		
		
	#show the dialog
	currSchema = Schema.GetCurrent()
	wgmlstExperType = bns.Database.ExperimentType(currSchema.WgMLSTExperTypeID)
	
	dlg = CalcStrainNamesDlg(wgmlstExperType, tmpNameFieldId, tmpNameHistoryFieldId, tmpCharViewId, tmpThresholds, tmpSplitNames, tmpMinPresenceThreshold)
	if not dlg.Show():
		bns.Stop()
		
	#store the current settings
	settings = ET.Element("hierarchical_strain_nomenclature_settings")
	ET.SubElement(settings, 'nameFieldId').text = dlg.selectedField
	ET.SubElement(settings, 'nameHistoryFieldId').text = dlg.selectedHistoryField
	ET.SubElement(settings, 'charViewId').text = dlg.selectedCharView
	ET.SubElement(settings, 'thresholds').text = '|'.join(str(f) for f in dlg.selectedThresholds)
	ET.SubElement(settings, 'splitNames').text = '1'if  dlg.selectedSplitNames else '0'
	ET.SubElement(settings, 'minPresenceThreshold').text = str(dlg.minPresenceThreshold)
	settings = bns.Database.Db.Info.SaveSetting(settingsKey, ET.tostring(settings), True)
		
	#calculate the damn thing
	allNamesFlName = os.path.join(bns.Database.Db.Info.SourceFilesDir, '9strain_naming_names_{}_{}.txt.gz'.format(wgmlstExperType.ID, dlg.selectedCharView))
	allProfilesFlName = os.path.join(bns.Database.Db.Info.SourceFilesDir, '9strain_naming_profiles_{}_{}.txt.gz'.format(wgmlstExperType.ID, dlg.selectedCharView))
	
	ok = os.path.isfile(allNamesFlName) and os.path.isfile(allProfilesFlName)
	if not ok:
		#put out a warning
		
		#delete files if only one is present
		if os.path.isfile(allNamesFlName): os.remove(allNamesFlName)
		if os.path.isfile(allProfilesFlName): os.remove(allProfilesFlName)
		
	calculator = Calculator(wgmlstExperType, dlg.selectedField, dlg.selectedHistoryField, dlg.selectedCharView, dlg.selectedThresholds, dlg.selectedSplitNames, dlg.minPresenceThreshold, allNamesFlName, allProfilesFlName)
	#calculator.DoCalc({})
	
	bns.Windows.BnsWindow(winId).StartAsyncCalculation(calculator.DoCalc, calculator.DoRefresh, async=False)
	


# inject section
try:
	bns.Windows.BnsWindow.AddCustomSection(WindowClass='main',  sectionID='WgsTools',  sectionName="WGS tools", parentSectionID='')
	bns.Windows.BnsWindow.AddCustomSection(WindowClass='main',  sectionID='CalcEngine',  sectionName="Calculation engine", parentSectionID='WgsTools')
	bns.Windows.BnsWindow.AddCustomSection(WindowClass='main',  sectionID='WgMlst',  sectionName="wgMLST", parentSectionID='WgsTools')
	bns.Windows.BnsWindow.AddCustomSection(WindowClass='main',  sectionID='WgMlstClient',  sectionName="wgMLST Client", parentSectionID='WgMlst')
	bns.Windows.BnsWindow.AddCustomSection(WindowClass='main',  sectionID='WgMlstCurator',  sectionName="wgMLST Curator", parentSectionID='WgMlst')
	bns.Windows.BnsWindow.AddCustomSection(WindowClass='main',  sectionID='WgSettings',  sectionName="Settings", parentSectionID='WgsTools')
except:
	pass

bns.Windows.BnsWindow.AddCustomSection(WindowClass='main',  sectionID='WgMlstStrainNaming', sectionName='wgMLST strain naming', parentSectionID='WgMlst')

# inject command
bns.Windows.BnsWindow.AddCustomCommand(
			WindowClass='main',
	   sectionID='WgMlstClient',
	   commandID='assign_strain_names',
	   caption='Assign strain names',
	   description='Assign strain names',
	   execute=Main
)

	
Main({})
	
	
