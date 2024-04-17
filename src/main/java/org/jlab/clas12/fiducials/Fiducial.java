package org.jlab.clas12.fiducials;

import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.JFrame;
import javax.swing.JTabbedPane;
import org.jlab.clas.physics.Particle;
import org.jlab.detector.base.DetectorType;
import org.jlab.detector.base.GeometryFactory;
import org.jlab.detector.geant4.v2.DCGeant4Factory;
import org.jlab.geom.base.ConstantProvider;
import org.jlab.geom.prim.Point3D;
import org.jlab.geom.prim.Trap3D;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.group.DataGroup;
import org.jlab.groot.base.GStyle;
import org.jlab.groot.data.DataLine;
import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.data.IDataSet;
import org.jlab.groot.data.TDirectory;
import org.jlab.groot.fitter.DataFitter;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.graphics.EmbeddedCanvasTabbed;
import org.jlab.groot.graphics.EmbeddedPad;
import org.jlab.groot.math.F1D;
import org.jlab.groot.math.Func1D;
import org.jlab.groot.math.FunctionFactory;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.*;
import org.jlab.logging.DefaultLogger;
import org.jlab.utils.benchmark.ProgressPrintout;
import org.jlab.utils.groups.IndexedList;
import org.jlab.utils.options.OptionParser;
        
        
/**
 *
 * @author devita
 */
public class Fiducial {
    // detector,  observable  
    Map<TabGroup, IndexedList<DataGroup>> dgs = new LinkedHashMap<>();
    
    double delta = 3;
    double[] thetas = {8, 15, 23};
    double dtheta = 0.5;
    int[] pids = {211, 2212};

    DCGeant4Factory dcFactory = null;
    
    Random rand = new Random();
    
    private static final Logger LOGGER = Logger.getLogger(Fiducial.class.getName());
    private static Level LEVEL = Level.CONFIG;
 
    
    public Fiducial() {
        this.createHistos();
        
        ConstantProvider dcProvider   = GeometryFactory.getConstants(DetectorType.DC, 11, "rga_fall2018");
        dcFactory = new DCGeant4Factory(dcProvider, DCGeant4Factory.MINISTAGGERON, false);

    }
    
    
    
    
    public final void createHistos() {

        GStyle.getAxisAttributesX().setTitleFontSize(24); //24
        GStyle.getAxisAttributesX().setLabelFontSize(18); //18
        GStyle.getAxisAttributesY().setTitleFontSize(24); //24
        GStyle.getAxisAttributesY().setLabelFontSize(18); //18
        GStyle.getAxisAttributesZ().setLabelFontSize(14); //14
        GStyle.setPalette("kDefault");
        GStyle.getAxisAttributesX().setLabelFontName("Arial");
        GStyle.getAxisAttributesY().setLabelFontName("Arial");
        GStyle.getAxisAttributesZ().setLabelFontName("Arial");
        GStyle.getAxisAttributesX().setTitleFontName("Arial");
        GStyle.getAxisAttributesY().setTitleFontName("Arial");
        GStyle.getAxisAttributesZ().setTitleFontName("Arial");
        GStyle.setGraphicsFrameLineWidth(1);
        GStyle.getH1FAttributes().setLineWidth(2);
        GStyle.getH1FAttributes().setOptStat("");

        DataGroup dgEB   = new DataGroup(2, DataType.values().length);
        for(DataType type : DataType.values()) {
            H2F hiptheta   = new H2F("hiptheta"+type.getName(), type.getName(), 200, 0, 6, 200, 5, 45);
            hiptheta.setTitleX("p (GeV)");
            hiptheta.setTitleY("#theta (deg)");
            H2F hiphitheta = new H2F("hiphitheta"+type.getName(), type.getName(), 200, -180, 180, 200, 5, 45);
            hiphitheta.setTitleX("#phi (deg)");
            hiphitheta.setTitleY("#theta (deg)");
            dgEB.addDataSet(hiptheta,   0+type.getId()*2);
            dgEB.addDataSet(hiphitheta, 1+type.getId()*2);
        }
        this.addDataGroup(TabGroup.KINEMATICS,  0, 0, 0, dgEB);
        for(int sector=0; sector<=6; sector++) {
            for(int layer=1; layer<=3; layer++) {
                for(int charge=1; charge>=-1; charge-=2) {
                    this.addDataGroup(TabGroup.DC, sector, layer, charge, this.createHistos(sector, layer, charge));
                }
            }
        }
    }
    
    private DataGroup createHistos(int sector, int layer, int charge) {
        DataGroup dg = new DataGroup(4, DataType.values().length);
        for(DataType type : DataType.values()) {
            double size = 30+layer*60;
            String title = type.getName()+"S"+sector+"L"+layer+"C"+charge;
            H2F hi1 = new H2F("hi"+type.getName()+"xy", title, 200, -size, size, 200, -size, size);  
            hi1.setTitleX("y (cm)"); 
            hi1.setTitleY("x (cm)"); 
            H2F hi2 = new H2F("hi"+type.getName()+"xycut", title, 200, -size, size, 200, -size, size);  
            hi2.setTitleX("y (cm)"); 
            hi2.setTitleY("x (cm)"); 
            H2F hi3 = new H2F("hi"+type.getName()+"thetaphi", title,200, -30, 30, 200, 0, 50);  
            hi3.setTitleX("#phi_t_r_a_j (deg)"); 
            hi3.setTitleY("#theta_t_r_a_j (deg)"); 
            H2F hi4 = new H2F("hi"+type.getName()+"thetaedge", title, 25, -10*0*layer, 10*layer, 30, 5, 35); 
    //                if(region==1)
    //                    hi5 = new H2F("hi"+type.getName()+region, "Region"+region, 40, -10, 10, 30, 5, 35);   
            hi4.setTitleX("edge (cm)"); 
            hi4.setTitleY("#theta_t_r_a_j (deg)"); 
            dg.addDataSet(hi1, 0+4*type.getId());
            dg.addDataSet(hi2, 1+4*type.getId());
            dg.addDataSet(hi3, 2+4*type.getId());
            dg.addDataSet(hi4, 3+4*type.getId());
        }
        return dg;
    }
    
    private void addDataGroup(TabGroup tab, int sector, int layer, int charge, DataGroup dg) {
        if(!dgs.containsKey(tab)) 
            dgs.put(tab, new IndexedList(3));
        
        IndexedList<DataGroup> list = dgs.get(tab);
        if(!list.hasItem(sector, layer))
            list.add(dg, sector, layer, charge);
    }
    
    public void processEvent(DataEvent event, int id, int selectedPid, int sec) {
        
        String type = DataType.getType(id).getName();
        
        DataBank recPart = null;
        DataBank recTraj = null;
        DataBank recTrac = null;
        DataBank mcPart = null;
        
        if (!event.hasBank("RUN::config")) {
            return;
        }
        int run = event.getBank("RUN::config").getInt("run", 0);

        if (event.hasBank("MC::Particle")) {
            mcPart = event.getBank("MC::Particle");
        }
        if (event.hasBank("REC::Particle")) {
            recPart = event.getBank("REC::Particle");
        }
        if (event.hasBank("REC::Traj")) {
            recTraj = event.getBank("REC::Traj");
        }
        if (event.hasBank("REC::Track")) {
            recTrac = event.getBank("REC::Track");
        }
        
        if (recPart != null && recTraj!=null && recTrac!=null) {
            if(recPart.getInt("pid", 0)!=11 || recPart.getShort("status", 0)>-2000) return;
            for(int i=0; i<recTraj.rows(); i++) {
                int pindex   = recTraj.getShort("pindex", i);
                int pid      = recPart.getInt("pid", pindex);
                int charge   = recPart.getByte("charge", pindex);
                int status   = recPart.getShort("status", pindex);
                int clas     = Math.abs(status)/1000;
                
                if(pindex==0 || clas!=2) continue;
                
                if(selectedPid==0 && pid==0) pid=charge*211;
                if(pid==0) continue;
                if(selectedPid!=0 && Math.abs(pid)!=selectedPid) continue;
                
                int detector = recTraj.getInt("detector", i);
                int layer    = recTraj.getInt("layer", i);
                double x     = recTraj.getFloat("x", i);
                double y     = recTraj.getFloat("y", i);
                double z     = recTraj.getFloat("z", i);
                double edge  = recTraj.getFloat("edge", i);
                int sector = 0;
                for(int j=0; j<recTrac.rows(); j++) {
                    if(recTrac.getShort("pindex", j)==pindex)
                        sector = recTrac.getByte("sector", j);
                }
                
                Particle part = new Particle(pid,
                                             recPart.getFloat("px", pindex),
                                             recPart.getFloat("py", pindex),
                                             recPart.getFloat("pz", pindex),
                                             recPart.getFloat("vx", pindex),
                                             recPart.getFloat("vy", pindex),
                                             recPart.getFloat("vz", pindex));
                if(part.p()<1 || charge==0) continue;
                
                Point3D trajTilted = new Point3D(x, y, z);
                Point3D trajLocal  = new Point3D(x, y, z);
                trajTilted.rotateZ(-Math.PI/3*(sector-1));
                trajLocal.rotateZ(Math.PI/2-Math.PI/3*(sector-1));
                if(layer==3 && detector==DetectorType.FTOF.getDetectorId())
                    trajLocal.rotateX(Math.toRadians(58.11));
                else
                    trajLocal.rotateX(Math.toRadians(25));
                
                if(detector==DetectorType.DC.getDetectorId() && (layer==6 || layer==18 || layer==36)) {
                    int region = (layer-1)/12+1;
                    double theta = Math.toDegrees(trajTilted.toVector3D().theta());
                    double phi   = Math.toDegrees(trajTilted.toVector3D().phi());
                                        
     //               System.out.println(region + " " + charge);
                    dgs.get(TabGroup.DC).getItem(sector,region, charge).getH2F("hi"+type+"xy").fill(trajLocal.x(), trajLocal.y());
                    if(edge>delta) dgs.get(TabGroup.DC).getItem(sector,region, charge).getH2F("hi"+type+"xycut").fill(trajLocal.x(), trajLocal.y());
                    dgs.get(TabGroup.DC).getItem(sector,region, charge).getH2F("hi"+type+"thetaphi").fill(phi, theta);
                    if(edge>0) {
                        int side = 1;//(int) Math.signum(trajLocal.x());
                        dgs.get(TabGroup.DC).getItem(sector,region, charge).getH2F("hi"+type+"thetaedge").fill(edge*side, theta);
                    }
                }
                else if(detector==DetectorType.HTCC.getDetectorId()) {
                    dgs.get(TabGroup.KINEMATICS).getItem(0,0,0).getH2F("hiptheta"+type).fill(part.p(), Math.toDegrees(part.theta()));
                    dgs.get(TabGroup.KINEMATICS).getItem(0,0, 0).getH2F("hiphitheta"+type).fill(Math.toDegrees(part.phi()), Math.toDegrees(part.theta()));
                }
            }		
        }

    }
    
    private void addHistos(DataGroup sum, DataGroup addendum) {
        int nrows = sum.getRows();
        int ncols = sum.getColumns();
        int nds   = nrows*ncols;
        for(int i=0; i<nds; i++) {
            List<IDataSet> dsList = sum.getData(i);
            for(IDataSet ds : dsList){
//                    System.out.println("\t --> " + ds.getName());
                if(ds instanceof H1F) {
                    H1F h1 = (H1F) ds;
                    H1F hi = addendum.getH1F(ds.getName());
                    h1.add(hi);
                }
                else if(ds instanceof H2F) {
                    H2F h1 = (H2F) ds;
                    H2F hi = addendum.getH2F(ds.getName());
                    h1.add(hi);
                }
            }
        }    
    }
    
    public void analyzeHisto(int nevents) {
        //create sector summed
        IndexedList<DataGroup> list = dgs.get(TabGroup.DC);
        for(int layer=1; layer<=3; layer++) {
            for(int charge=1; charge>=-1; charge-=2) {
                if(!list.hasItem(0, layer, charge))
                    continue;
                DataGroup dg0 = list.getItem(0, layer, charge);
                for(int is=0; is<6; is++) {
                    DataGroup dgi = list.getItem(is+1, layer, charge);
                    this.addHistos(dg0,dgi);
                }
            }
        }
        double[] edge0 = {3, 1, 8};
        for(int sector=0; sector<=6; sector++) {
            for(int layer=1; layer<=3; layer++) {
                for(int charge=1; charge>=-1; charge-=2) {
                    for(int itype=0; itype<DataType.values().length; itype++) {
                        String type = DataType.getType(itype).getName();
                        int icol = (itype+1)*2;
                        H2F h2 = dgs.get(TabGroup.DC).getItem(sector,layer,charge).getH2F("hi"+type+"thetaedge");
                        dgs.get(TabGroup.DC).add(this.fitEdges(h2, icol, edge0[layer-1]), sector, layer+10*(itype+1), charge);
                    }
                }
            }
        }
    }
    
    private DataGroup fitEdges(H2F h2, int icol, double edge0) {
        DataGroup dg = new DataGroup(6, 5);
        List<H1F> slices = h2.getSlicesY();
        for(int itheta=0; itheta<slices.size(); itheta++) {
            H1F h1 = slices.get(itheta);
            h1.setLineColor(icol);
            h1.setName(h2.getName()+"_"+h2.getDataY(itheta));
            h1.setTitle("#theta="+h2.getDataY(itheta)+"(deg)");
            double xmax = h1.getDataX(h1.getDataSize(0)-1)*0.8;
            double dx = h1.getDataX(1)-h1.getDataX(0);
            edge0 = Fiducial.getThresholdCrossing(h1, 0.8);
            if(Fiducial.getIntegralIDataSet(h1, xmax-dx, xmax+dx)<5) continue;
            if(Fiducial.getIntegralIDataSet(h1, xmax, xmax+2*dx) < Fiducial.getIntegralIDataSet(h1, xmax-2*dx, xmax)*0.8) continue;
            EdgeSlice f1 = new EdgeSlice("f1"+h2.getName()+"_"+h2.getDataY(itheta), 0, xmax);
            f1.setParameter(0, h1.getMax());
            f1.setParameter(1, edge0);
            f1.setParameter(2, 0.5);
            f1.setParameter(3, 0);
            f1.setParLimits(2,   0.0,  5.0);
            f1.setParLimits(3, -20.0, 20.0);
            f1.setLineWidth(2);
            DataFitter.fit(f1, h1, "Q");
            if(f1.isFitValid() || true) {
                dg.addDataSet(h1, itheta); 
            }
        }
        return dg;
    }
 
        
    public class EdgeSlice extends Func1D{

            public EdgeSlice(String name, double min, double max) {
                    super(name, min, max);
                    super.addParameter("amp1");
                    super.addParameter("mean1");
                    super.addParameter("sigma1");
                    super.addParameter("slope1");
//                    super.addParameter("amp2");
//                    super.addParameter("mean2");
//                    super.addParameter("sigma2");
//                    super.addParameter("slope2");
            }

            @Override
            public double evaluate(double x){
                    double value = 0.0;
                    
//                    if(x>0) {
                        if(Math.abs(x)<this.getParameter(1)) 
                            value = this.getParameter(0)*FunctionFactory.gauss(Math.abs(x),this.getParameter(1), this.getParameter(2));
                        else
                            value = this.getParameter(0)+(Math.abs(x)-this.getParameter(1))*this.getParameter(3);
//                    }
//                    else {
//                        if(x>this.getParameter(5)) 
//                            value = this.getParameter(4)*FunctionFactory.gauss(x,this.getParameter(5), this.getParameter(6));
//                        else
//                            value = this.getParameter(4)+(x-this.getParameter(5))*this.getParameter(7);
//                    }
                    return value;
            }
    }


    private DataGroup readDataGroup(String folder, TDirectory dir, DataGroup dg) {
        int nrows = dg.getRows();
        int ncols = dg.getColumns();
        int nds   = nrows*ncols;
        DataGroup newGroup = new DataGroup(ncols,nrows);
        for(int i = 0; i < nds; i++){
            List<IDataSet> dsList = dg.getData(i);
            for(IDataSet ds : dsList){
                if(dir.getObject(folder, ds.getName())!=null) {
                    IDataSet dsread = dir.getObject(folder, ds.getName());
                    if(dsread instanceof H1F) {
                        H1F h1 = (H1F) dsread;
                        Func1D f1 = (Func1D) dir.getObject(folder, "f"+h1.getName());
                        if(f1!=null)
                            h1.setFunction(f1);
                    }
                    if(dsread instanceof H1F && ((H1F) dsread).getFunction()!=null) {
                        Func1D dsf = ((H1F) dsread).getFunction();
                        dsf.setLineColor(2);
                        dsf.setLineWidth(2);
                        String opts = "11";
                        for(int k=0; k<dsf.getNPars(); k++) opts += "1";
                        dsf.setOptStat(opts);
                    }
                    else if(dsread instanceof F1D && dg.getF1D(dsread.getName())!=null) {
                        F1D dgf = (F1D) dg.getF1D(dsread.getName());
                        F1D dsf = (F1D) dsread;
                        dsf.setLineColor(dgf.getLineColor());
                        dsf.setLineWidth(dgf.getLineWidth());
                        dsf.setOptStat(dgf.getOptStat());
                    }
                    newGroup.addDataSet(dsread,i);                        
                }
            }
        }
        return newGroup;
    }

    public void readHistos(String fileName, String optStats) {
        LOGGER.log(LEVEL,"Opening file: " + fileName);
        TDirectory dir = new TDirectory();
        dir.readFile(fileName);
        for(TabGroup tab : dgs.keySet()) {
            String folder = tab.getName();
            IndexedList<DataGroup> list = dgs.get(tab);
            for(int sector=0; sector<=6; sector++) {
                for(int layer=0; layer<=3; layer++) {
                    for(int charge=1; charge>=-1; charge-=2) {
                        if(!list.hasItem(sector, layer, charge))
                            continue;
                        String name = "S" + sector + "L" + layer + "C" + charge;
                        dir.cd("/" + folder + "/" + name);
                        list.add(this.readDataGroup(folder + "/" + name, dir, list.getItem(sector, layer, charge)), sector, layer, charge);
                    }
                }
            }
        }
    }

    public void saveHistos(String fileName) {
        LOGGER.log(LEVEL,"\nSaving histograms to file " + fileName);
        TDirectory dir = new TDirectory();
        for(TabGroup tab : dgs.keySet()) {
            String folder = tab.getName();
            dir.cd();
            dir.mkdir("/" + folder);
            dir.cd("/" + folder);
            IndexedList<DataGroup> list = dgs.get(tab);
            for(int sector=0; sector<=6; sector++) {
                for(int layer=0; layer<=3; layer++) {
                    for(int charge=1; charge>=-1; charge-=2) {
                        if(!list.hasItem(sector, layer, charge))
                            continue;
                        String name = "S" + sector + "L" + layer + "C" + charge;
                        dir.mkdir("/" + folder + "/" + name);
                        dir.cd("/" + folder + "/" + name);
                        System.out.println("Saving histograms in " + "/" + folder + "/" + name);
                        this.writeDataGroup(dir, list.getItem(sector, layer));
                    }
                }
            }
        }
        dir.writeFile(fileName);
    }        

    private void writeDataGroup(TDirectory dir, DataGroup dg) {
        int nrows = dg.getRows();
        int ncols = dg.getColumns();
        int nds   = nrows*ncols;
        for(int i=0; i<nds; i++) {
            List<IDataSet> dsList = dg.getData(i);
            for(IDataSet ds : dsList){
//                    System.out.println("\t --> " + ds.getName());
                dir.addDataSet(ds);
                if(ds instanceof H1F) {
                    H1F h1 = (H1F) ds;
                    if(h1.getFunction()!=null) {
                        Func1D f1 = h1.getFunction();
                        f1.setName("f"+h1.getName());
                        dir.addDataSet(f1);
                    }
                }
            }
        }
    }
        
    private static double getEfficiency(H1F hi) {
        
        double phiRef = 2;
        int nRef = 500;
        
        if(Fiducial.getIntegralIDataSet(hi, -phiRef, phiRef)<nRef)
            return 0;
        else {
            double average = Fiducial.getAverageIDataSet(hi, -phiRef, phiRef);
            double phiMin = Fiducial.getThresholdCrossing(hi, 0.3, 1);
            double phiMax = Fiducial.getThresholdCrossing(hi, 0.3,-1);
            double phiThr = (phiMax-phiMin)/2;
            double integral = Fiducial.getIntegralIDataSet(hi, -phiThr, phiThr);
            double efficiency = integral/average/phiThr/2;
            System.out.println(average + " " + phiMin + " " + phiMax + " " + integral + " " + efficiency);
            if(average<nRef/10)
                return 0;
            else if(phiThr<phiRef*2)
                return 0;
            else
                return efficiency;
        }
    }
        
    public static int getMaximumBinBetween(H1F histo, double min, double max) { 
        int nbin = histo.getData().length;
        double x_val_temp;
        double x_val;
        double y_max_temp;
        double y_max = 0;
        int max_bin_num = histo.getMaximumBin();
        for (int i = 0; i < nbin; i++) { 
            x_val_temp = histo.getAxis().getBinCenter(i);
            if (x_val_temp >= min && x_val_temp <= max) {
                y_max_temp = histo.getBinContent(i);
                if (y_max_temp > y_max) {
                    y_max = y_max_temp;
                    max_bin_num = i;
                }
            }
        }
        return max_bin_num;
    }
    
    private static double getIntegralIDataSet(IDataSet data, double min, double max) {
            double nEntries = 0;
            for (int i = 0; i < data.getDataSize(0); i++) {
                    double x = data.getDataX(i);
                    double y = data.getDataY(i);
                    if (x > min && x < max && y != 0) {
                            nEntries += y;
                    }
            }
            return (double) nEntries*(data.getDataX(1)-data.getDataX(0));
    }

    private static double getAverageIDataSet(IDataSet data, double min, double max) {
            int nsamples = 0;
            double sum = 0;
            for (int i = 0; i < data.getDataSize(0); i++) {
                    double x = data.getDataX(i);
                    double y = data.getDataY(i);
                    if (x > min && x < max && y != 0) {
                            nsamples++;
                            sum += y;
                    }
            }
            return sum / nsamples;
    }

    private static double getMeanIDataSet(IDataSet data, double min, double max) {
            int nsamples = 0;
            double sum = 0;
            double nEntries = 0;
            for (int i = 0; i < data.getDataSize(0); i++) {
                    double x = data.getDataX(i);
                    double y = data.getDataY(i);
                    if (x > min && x < max && y != 0) {
                            nsamples++;
                            sum += x * y;
                            nEntries += y;
                    }
            }
            return sum / (double) nEntries;
    }
   
    private static GraphErrors getThresholdCrossingProfile(H2F hi, double fraction) {
        GraphErrors graph = new GraphErrors();
        graph.setTitleX(hi.getTitleX());
        graph.setTitleY(hi.getTitleY());
        List<H1F> hslices = hi.getSlicesX();
        for(int ix=0; ix<hi.getDataSize(0); ix++) {
            H1F hix = hslices.get(ix);
            if(hix.getMax()>10) {
                graph.addPoint(hi.getDataX(ix), Fiducial.getThresholdCrossing(hix, fraction), 0, (hix.getDataX(1)-hix.getDataX(0))/2);
            }
        }
        return graph;
    }
    
    private static double getThresholdCrossing(H1F hi, double fraction) {
        return Fiducial.getThresholdCrossing(hi, fraction, 1);
    }
    
    private static double getThresholdCrossing(H1F hi, double fraction, int dir) {
        if(fraction<0 || fraction>=1) {
            System.out.println("[ERROR] Invalid constant fraction threshold " + fraction + " for histogram " + hi.getName());
            System.exit(1);
        }
        double hiMax = hi.getMax();
        int iHalf = -1;
        for(int i=0; i<hi.getDataSize(0); i++) {
            int bin = dir>0 ? i : hi.getDataSize(0)-1-i;
            if(hi.getDataY(bin)>hiMax*fraction) {
                iHalf = bin;
                break;
            }
        }
        return hi.getDataX(iHalf);
    }
    private static double getRMSIDataSet(IDataSet data, double min, double max) {
            int nsamples = 0;
            double mean = getMeanIDataSet(data, min, max);
            double sum = 0;
            double nEntries = 0;

            for (int i = 0; i < data.getDataSize(0); i++) {
                    double x = data.getDataX(i);
                    double y = data.getDataY(i);
                    if (x > min && x < max && y != 0) {
                            nsamples++;
                            sum += Math.pow(x - mean, 2) * y;
                            nEntries += y;
                    }
            }
            return Math.sqrt(sum / (double) nEntries);
    }
        
    private EmbeddedCanvasTabbed getGroupByLayer(TabGroup tab, int order) {
        EmbeddedCanvasTabbed canvas = null;        
        for(int sector=0; sector<=6; sector++) {
            DataGroup dn = new DataGroup(3, 2*DataType.values().length);
            String name = sector==0 ? "All sectors" : "Sector " + sector;
            if(canvas==null) {
                canvas = new EmbeddedCanvasTabbed(name);
            }
            else {
                canvas.addCanvas(name);
            }
            for(int layer=1; layer<=3; layer++) {
                for(int charge=1; charge>=-1; charge-=2) {
                    DataGroup dg = dgs.get(tab).getItem(sector, layer, charge);
                    for(int itype=0; itype<DataType.values().length; itype++) {
                        List<IDataSet> ds = dg.getData(order+itype*dg.getColumns());
                        for(IDataSet d : ds) {
                            dn.addDataSet(d, layer-1+3*((1-charge)/2+2*itype));
                        }
                    }
                }
            }
            canvas.getCanvas(name).draw(dn);
            if(order<2) this.drawDC(canvas.getCanvas(name));
            for(EmbeddedPad pad : canvas.getCanvas(name).getCanvasPads()) {
                pad.setTitleFont("Arial");
                pad.setTitleFontSize(18);
                if(!pad.getDatasetPlotters().isEmpty() && pad.getDatasetPlotters().get(0).getDataSet() instanceof H2F) 
                    pad.getAxisZ().setLog(true);
                if(!pad.getDatasetPlotters().isEmpty() && pad.getDatasetPlotters().get(0).getDataSet() instanceof GraphErrors) {
                    pad.getAxisX().setRange(5, 36);
                    pad.getAxisY().setRange(0, 10);
                }
            }
        }
        return canvas;
    }
    
    private EmbeddedCanvasTabbed getEdgeGraphs(TabGroup tab) {
        EmbeddedCanvasTabbed canvas = new EmbeddedCanvasTabbed("Edges");
        DataGroup dn = new DataGroup(3,2);
        for(int sector=0; sector<=6; sector++) {
            for(int layer=1; layer<=3; layer++) {
                for(int charge=1; charge>=-1; charge-=2) {
                    for(int itype=0; itype<DataType.values().length; itype++) {
                        DataGroup dg = dgs.get(tab).getItem(sector, layer+10*(itype+1), charge);
                        int icol = 2+2*itype;
                        String name = DataType.getType(itype).getName()+"S"+sector+"L"+layer+"C"+charge;
                        GraphErrors mean = new GraphErrors("EdgeMean"+name);
                        mean.setTitleX("#theta (deg)");
                        mean.setTitleY("edge (cm)");
                        mean.setMarkerStyle((1-charge)/2);
                        mean.setMarkerSize(5);
                        mean.setMarkerColor(icol+(1-charge)/2*5);
                        GraphErrors sigma = new GraphErrors("EdgeSigma"+name);
                        sigma.setTitleX("#theta (deg)");
                        sigma.setTitleY("#sigma(edge) (cm)");
                        sigma.setMarkerStyle((1-charge)/2);
                        sigma.setMarkerSize(5);
                        sigma.setMarkerColor(icol+(1-charge)/2*5);
                        for(int i=0; i<dg.getColumns()*dg.getRows(); i++) {
                            List<IDataSet> ds = dg.getData(i);
                            for(IDataSet d : ds) {
                                if(d instanceof H1F && ((H1F) d).getFunction()!=null) {
                                    Func1D f1 = ((H1F) d).getFunction();
                                    String theta = f1.getName().split("_")[1];
                                    mean.addPoint(Double.parseDouble(theta), f1.getParameter(1), 0, f1.parameter(1).error());
                                    sigma.addPoint(Double.parseDouble(theta), f1.getParameter(2), 0, f1.parameter(2).error());
                                }
                            }
                        }
                        if(mean.getDataSize(0)>0) {
                            dn.addDataSet(mean,  layer-1);
                            dn.addDataSet(sigma, layer+2);
                        }
                    }
                }
            }
        }
        canvas.getCanvas().draw(dn);
        return canvas;
    }
    
    private EmbeddedCanvasTabbed getEdgeGroup(TabGroup tab) {
        EmbeddedCanvasTabbed canvas = null;
        for(int sector=0; sector<=6; sector++) {
            for(int layer=1; layer<=3; layer++) {
                String name = sector==0 ? "All sectors" : "Sector " + sector;
                name += " Layer " + layer;
                if(canvas==null) {
                    canvas = new EmbeddedCanvasTabbed(name);
                }
                else {
                    canvas.addCanvas(name);
                }
                for(int charge=1; charge>=-1; charge-=2) {
                    canvas.getCanvas().draw(dgs.get(tab).getItem(sector, layer+10, charge));
                    canvas.getCanvas().draw(dgs.get(tab).getItem(sector, layer+20, charge));
                }
            }
        }
        return canvas;
    }
    
    private JTabbedPane getCanvas() {
        JTabbedPane panel = new JTabbedPane();
        for(TabGroup tab : dgs.keySet()) {
            IndexedList<DataGroup> list = dgs.get(tab);
            if(tab==TabGroup.KINEMATICS) {
                EmbeddedCanvasTabbed canvas = new EmbeddedCanvasTabbed(tab.getName());
                canvas.getCanvas().draw(list.getItem(0, 0, 0));
                panel.add(tab.getName(), canvas);
                for(EmbeddedPad pad : canvas.getCanvas().getCanvasPads()) {
                    pad.setTitleFont("Arial");
                    pad.setTitleFontSize(18);
                    if(!pad.getDatasetPlotters().isEmpty() && pad.getDatasetPlotters().get(0).getDataSet() instanceof H2F) 
                        pad.getAxisZ().setLog(true);
                }        
            }
            else {
                for(int order=0; order<6; order++) {
                    String name = tab.getName()+order;
                    if(order<4)
                        panel.add(name, this.getGroupByLayer(tab, order));
                    else if(order==4)
                        panel.add(name, this.getEdgeGroup(tab));
                    else
                        panel.add(name, this.getEdgeGraphs(tab));
                }   
            }
        }
        return panel;
    }
    
    private void drawDC(EmbeddedCanvas canvas) {
        int nrow = canvas.getNRows();
        for(int layer : new int[]{6, 18, 36}) {
            int superlayer = layer/6;
            for(int irow=0; irow<nrow; irow++) {
                canvas.cd((superlayer-1)/2+irow*3);
                Trap3D contour = dcFactory.getTrajectorySurface(0, superlayer-1, (layer-1)%6);
                contour.rotateZ(Math.PI/2);
                contour.rotateX(Math.toRadians(25));
                DataLine[] lines = new DataLine[4];
                for(int i=0; i<lines.length; i++) {
                    int i1 = i;
                    int i2 = i+1; if(i2==lines.length) i2=0;
                    lines[i] = new DataLine(contour.point(i1).x(),contour.point(i1).y(),contour.point(i2).x(),contour.point(i2).y());
                    lines[i].setLineColor(2);
                    canvas.draw(lines[i]);
                }
            }
        }
    }

    public static enum TabGroup {
        
        UNDEFINED  ( 0, "Undefined"),
        KINEMATICS ( 1, "Kinematics"),    
        DC         ( 2, "DC");
        
        private final int tabId;
        private final String tabName;
        
        TabGroup(){
            tabId = 0;
            tabName = "UNDEFINED";
        }
        
        TabGroup(int id, String name) {
            tabId = id;
            tabName = name;
        }
    
        public String getName() {
            return tabName;
        }

        public int getId() {
            return tabId;
        }

        public static TabGroup getType(String name) {
            name = name.trim();
            for(TabGroup tab: TabGroup.values())
                if (tab.getName().equalsIgnoreCase(name)) 
                    return tab;
            return UNDEFINED;
        }
        
        public static TabGroup getType(Integer id) {

            for(TabGroup tab: TabGroup.values())
                if (tab.getId() == id) 
                    return tab;
            return UNDEFINED;
        }
    }
    
    public static enum DataType {
        
        DATA      ( 0, "Data"),    
        MC        ( 1, "MC");
        
        private final int typeId;
        private final String typeName;
        
        DataType(){
            typeId = 0;
            typeName = "UNDEFINED";
        }
        
        DataType(int id, String name) {
            typeId = id;
            typeName = name;
        }

        public String getName() {
            return typeName;
        }

        public int getId() {
            return typeId;
        }

        public static DataType getType(String name) {
            name = name.trim();
            for(DataType type: DataType.values())
                if (type.getName().equalsIgnoreCase(name)) 
                    return type;
            return null;
        }
        
        public static DataType getType(Integer id) {
            for(DataType type: DataType.values())
                if (type.getId() == id) 
                    return type;
            return null;
        }
    }
    
    public static void main(String[] args) {
        

        OptionParser parser = new OptionParser("fiducials");
        parser.setRequiresInputList(true);
        parser.addOption("-n","-1", "maximum number of events to process");
        parser.addOption("-o","", "output file name");
        parser.addOption("-p","0", "particlr id (0=all)");
        parser.addOption("-s","0", "sector id (0=all)");
        parser.addOption("-d", "1", "open graphical window (1) or run in batch mode (0)");
        parser.addOption("-r", "0", "read histogram file");
        parser.parse(args);
        
        
        int     maxEvents = parser.getOption("-n").intValue();
        int     pid       = parser.getOption("-p").intValue();
        int     sector    = parser.getOption("-s").intValue();
        boolean display   = parser.getOption("-d").intValue()==1;
        boolean read      = parser.getOption("-r").intValue()==1;
        String  histos    = parser.getOption("-o").stringValue();
        
        if(!display) System.setProperty("java.awt.headless", "true");
        DefaultLogger.debug();
        
        Fiducial analysis = new Fiducial();        

        List<String> inputFiles = parser.getInputList();
        
        ProgressPrintout progress = new ProgressPrintout();

        if(!read) {
            for(int type=0; type<inputFiles.size(); type++) {
                int counter=-1;
                String input = inputFiles.get(type);
                HipoDataSource  reader = new HipoDataSource();
                reader.open(input);
                while(reader.hasEvent()) {

                    counter++;

                    DataEvent ev = reader.getNextEvent();

                    analysis.processEvent(ev, type, pid, sector);

                    progress.updateStatus();

                    if(maxEvents>0){
                        if(counter>=maxEvents) break;
                    }

                }
                progress.showStatus();
                reader.close();
            }   
        }
        else {
            analysis.readHistos(inputFiles.get(0), "");
        }
        analysis.analyzeHisto(0);
        if(!histos.isBlank())
            analysis.saveHistos(histos);
        if(display) {
            JFrame frame = new JFrame("Fiducials");
            frame.setSize(1500,1000);
            frame.add(analysis.getCanvas());
            frame.setLocationRelativeTo(null);
            frame.setVisible(true);     
        }
    }
}



