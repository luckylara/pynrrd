import numpy as np
from collections import OrderedDict
import os.path as path
import gzip

# Author: David Qixiang Chen
# email: qixiang.chen@gmail.com
#
# utility nrrd header reader for .nhdr
class NrrdHeader:
    def fromNiftiHeader(nifti_header):
        header = NrrdHeader()
        aff = nifti_header.get_best_affine()

        header.correctSpaceRas()
        header.setValue('space directions', aff[:3,:3])
        header.setValue('space origin', aff[:3,3])
        header.setValue('sizes', nifti_header.get_data_shape())
        header.setValue('dimension', str(len(nifti_header.get_data_shape())))

        return header

    fromNiftiHeader = staticmethod(fromNiftiHeader)


    b0num = 0
    def __init__ (self, dic=None):        
        # default headers
        self._data = OrderedDict({
            "header": "NRRD0004", 
            "type": "double", 
            "dimension": "3", 
            "space": "left-posterior-superior", 
            "sizes": [0,0,0], 
            "space directions": [[1, 0.0, 0.0],[0.0, -1, 0.0], [0.0, 0.0, 1]], 
            "kinds": ["domain", "domain", "domain"], 
            "endian": "little", 
            "encoding": "gzip", 
            "space origin": [0, 0, 0]
        })
           
        if dic:
            self._data = dic
        if dic.has_key('b0num'):
            self.b0num = self._data['b0num']

    def getAffine(self):
        affine = np.eye(4)
        affine[:3,:3] = np.sign(self['space directions'][:3])
        affine[:3,3] = self['space origin']

        return affine


    def getDwiGradients(self):
        return self._data['DWMRI_gradient']

    def setDwiGradients(self, vec):
         self._data['DWMRI_gradient'] = vec

    def setValue(self,key, val):
        self._data[key] = val

    def __getitem__(self, key):
        if key in self._data:          
            return self._data[key]
        return False

    def __setitem__(self, key, val):
        self._data[key] = val

    def getKeys(self):
        return self._data.keys()

    def isDTMR(self):
        return self.b0num > 0

    def getBval(self):
        return float(self['DWMRI_b-value'])

    def getBvals(self):
        b0 = self.getBval()
        bvals = np.repeat(float(b0), len(self.getDwiGradients()))
        return bvals

    def correctSpaceRas(self):
        def print_info(self):
            spaceraw = self['space']
            space = spaceraw.split('-')
            print 'space:',space
            origin = self['space origin']
            print 'origin:',origin
            directions = self['space directions']
            print 'space directions:',directions
            frame = np.array(self['measurement frame'])
            print 'measurement frame:',frame
        # don't invert dwi vectors, slicer doesn't do this internally
        #dwivec = self.getDwiGradients()
        print '==============before============'
        print_info(self)

        spaceraw = self['space']
        space = spaceraw.split('-')
        origin = self['space origin']
        directions = self['space directions']
        frame = np.array(self['measurement frame'])

        nospace = directions.index(np.nan)
        directions.remove(np.nan)
        dvec = np.array([directions])

        if space[0] == 'left':
            print '---------------------------'
            print 'invert left to right'
            space[0] = 'right'
            origin[0] *= -1
            dvec[:,0] *= -1
            frame[:,0] *= -1
            #dwivec[:,0] *= -1

        if space[1] == 'posterior':
            print '---------------------------'
            print 'invert posterior to anterior'
            space[1] = 'anterior'
            origin[1] *= -1
            dvec[:,1] *= -1
            frame[:,1] *= -1
            #dwivec[:,1] *= -1

        if space[1] == 'inferior':
            print '---------------------------'
            print 'invert inferior to superior'
            space[2] = 'superior'
            origin[2] *= -1
            dvec[:,2] *= -1
            frame[:,2] *= -1
            #dwivec[:,2] *= -1

        self.setValue('space','-'.join(space))
        dvec = np.insert(dvec.astype(np.object),nospace,np.nan,1).tolist()[0]
        dvec[nospace] = np.nan
        self.setValue('space directions',dvec)
        self.setValue('measurement frame', frame.tolist())
        #self.setDwiGradients(dwivec)

        print '==============after============'
        print_info(self)
        #dwivec = self.getDwiGradients()
        #print dwivec



class NrrdReader:
    grdkey = 'DWMRI_gradient'
    b0num = 'b0num'

    def load(self, filename, get_raw=False):
        filedir = path.dirname(path.abspath(filename))

        TFILE = open(filename, 'r')
        strbuf =""
        bindata = None
        if filename.find("nrrd") > -1:
            data = TFILE.read().split("\n\n",1)
            strbuf = data[0].split('\n')
            bindata = data[1]

        else:
            strbuf = TFILE.read().split('\n')

        TFILE.close()



        params = OrderedDict()

        for line in strbuf:
            if len(line) == 0:
                break;
            line = line.strip()
            if line == '':
                continue

            if line.startswith('#'):
                continue

            if line.startswith('NRRD'):
                params['header'] = line
            else:
                key,val = line.split(':')
                key = key.strip()
                if key=='space' and params.has_key('space') and params['space']!='':
                    continue

                val = val.replace('=','').strip()


                if key == 'sizes':
                    val = self.getVals(val, 'int')
                elif key.startswith(self.grdkey):
                    val = self.getVals(val, 'float')
                elif key.startswith('space origin'):
                    val = self.getVals(val, 'float')
                elif key.startswith('space directions'):
                    val = self.getVals(val, 'float')
                elif key.startswith('measurement frame'):
                    val = self.getVals(val, 'float')
                elif key.startswith('data file'):
                    val = val.strip()
                else:
                    val = self.getVals(val)

                if key.startswith(self.grdkey):
                    if not params.has_key(self.grdkey):
                        params[self.grdkey] = []
                    if not params.has_key(self.b0num):
                        params[self.b0num] = 0

                    params[self.grdkey].append(val)

                    if val[0] == 0 and val[1] == 0 and val[2] == 0:
                        #only want one b0 line, very hacky
                        #if params[self.b0num] < 1:
                        #    params[self.grdkey].append(val)                        
                        params[self.b0num] += 1
                        
                    #else:
                    #    params[self.grdkey].append(val)

                else:
                    params[key]= val

        if (params.has_key(self.grdkey)):
            dwivec = params[self.grdkey]
            params[self.grdkey] = dwivec

        params = NrrdHeader(params)
        if not get_raw:
            if not bindata and params['data file']:
                bin_file = params['data file']
                r_func = open
                if params['encoding'] == 'gzip':
                    r_func = gzip.open

                with r_func(path.join(filedir, bin_file), 'rb') as fp:
                    bindata = fp.read()
            
            elif bindata and params['encoding'] == 'gzip':
                from StringIO import StringIO
                bindata = gzip.GzipFile(fileobj=StringIO(bindata)).read()

            if bindata:
                type_val = params['type']
                if  type_val == 'short':
                    m_type = np.short
                elif type_val == 'float':
                    m_type = np.float32
                elif type_val == 'double':
                    m_type = np.double

                a = np.fromstring(bindata, dtype=m_type)

                bindata = np.reshape(a, params['sizes'])
        return params, bindata

    def asDtype(self, dtype, value):
        if dtype=='float':
            try:
                val = float(value)
            except ValueError:
                val = np.nan
            return val
        elif dtype=='int':
            return int(value)
        else:
            return value

    def getVals(self,input, dtype='str'):
        input = input.strip()
        if input.find('(') > -1:
            input = input.replace('(','|')
            input = input.replace(')','|')
            input = input.split('|')
        else:
            input = input.split(' ')

        res = []
        for i in input:
            row=[]
            i = i.strip()
            if i == '':
                continue
            if i.find(',') > -1:
                i = i.split(',')
                for k in range(len(i)):
                    row.append(self.asDtype(dtype, i[k].strip()))
                res.append(row)
            else:
                row = self.asDtype(dtype, i)
                res.append(row)
        if len(res) == 1:
            a = [res]
            res = res[0]
        return res

class NrrdWriter:
    def write(self,nrrdheader, output):
        FILE = open(output,'w')

        keys = nrrdheader.getKeys()
        eq = ': '
        for k in keys:
            is_write=True
            val = nrrdheader[k]
            if k=='header':
                line=val+'\n'
            elif k=='b0num':
                continue
            elif k=='modality' or k=='DWMRI_b-value':
                eq = ':='
                line = "%s%s%s\n" % (k,eq,val)
            elif k=='DWMRI_gradient':
                c = 0
                for i in val:
                    i = self.formatOutput(i, brac=False, dim=' ')
                    line = "%s_%04d:=%s\n" % (k, c, i)
                    c+=1
                    print line,
                    FILE.write(line)
                    is_write=False
            else:
                if k=='space origin':
                    val = self.formatOutput(val)
                elif k=='space directions':
                    val = self.formatOutput(val, nested=True)
                elif k=='measurement frame':
                    val = self.formatOutput(val, nested=True)
                else:
                    val = self.formatOutput(val,brac=False, dim=' ')
                line = "%s%s%s\n" % (k,eq,val)
            if is_write:
                print line,
                FILE.write(line)

        FILE.close()
    def formatOutput(self, lis, nested=False, brac=True, dim=','):
        res=''
        if type(lis) is str:
            return lis

        if not nested:
            if brac:
                res='('

            for i in lis:
                res+=str(i)+dim

            res = res[:-1]

            if brac:
                res+=')'
        else:
            for i in lis:
                if i is np.nan:
                    t=' none '
                else:
                    t='('
                    for j in i:
                        t+=str(j)+dim
                    t=t[:-1]
                    t+=') '
                res += t
        return res



