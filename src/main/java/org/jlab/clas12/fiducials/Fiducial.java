package org.jlab.clas12.fiducials;

import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import javax.swing.JFrame;
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
import org.jlab.utils.options.OptionParser;
        
        
/**
 *
 * @author devita
 */
public class Fiducial {
    
    Map<TabGroup, Map<Integer,DataGroup>> dgs = new LinkedHashMap<>();
    
    double delta = 3;
    double[] thetas = {8, 15, 23};
    double dtheta = 0.5;
    int[] pids = {211, 2212};

    DCGeant4Factory dcFactory = null;
    
    Random rand = new Random();
    
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
        DataGroup dgTraj = new DataGroup(3, 2);
        DataGroup dgPhi  = new DataGroup(3, 2);
        DataGroup dgCut  = new DataGroup(3, 2);
        DataGroup dgEdge = new DataGroup(3, 2);
        DataGroup dgPhiT = new DataGroup(3, 2);
        DataGroup dgEdg1 = new DataGroup(8, 6);
        DataGroup dgEdg2 = new DataGroup(8, 6);
        DataGroup dgEdg3 = new DataGroup(8, 6);
        DataGroup dgSecT = new DataGroup(3, 2);
        DataGroup dgEff  = new DataGroup(3, 2);
        for(DataType type : DataType.values()) {
            H2F hiptheta   = new H2F("hiptheta"+type.getName(), type.getName(), 200, 0, 6, 200, 5, 45);
            hiptheta.setTitleX("p (GeV)");
            hiptheta.setTitleY("#theta (deg)");
            H2F hiphitheta = new H2F("hiphitheta"+type.getName(), type.getName(), 200, -180, 180, 200, 5, 45);
            hiphitheta.setTitleX("#phi (deg)");
            hiphitheta.setTitleY("#theta (deg)");
            dgEB.addDataSet(hiptheta,   0+type.getId()*2);
            dgEB.addDataSet(hiphitheta, 1+type.getId()*2);
            for(int i=0; i<3; i++) {
                int region = i+1;
                double size = 90+i*60;
                H2F hi1 = new H2F("hi"+type.getName()+region, "Region "+region, 200, -size, size, 200, -size, size);  
                hi1.setTitleX("y (cm)"); 
                hi1.setTitleY("x (cm)"); 
                H2F hi2 = new H2F("hi"+type.getName()+region, "Region"+region, 200, -30, 30, 200, 0, 50);  
                hi2.setTitleX("#phi_t_r_a_j (deg)"); 
                hi2.setTitleY("#theta_t_r_a_j (deg)"); 
                H2F hi3 = new H2F("hi"+type.getName()+region, "Region"+region, 200, -size, size, 200, -size, size);  
                hi3.setTitleX("y (cm)"); 
                hi3.setTitleY("x (cm)"); 
                for(int k=0; k<4; k++) {
                    H1F hi4 = new H1F("hi"+type.getName()+region+k, "Region"+region, 100, -30, 30);  
                    hi4.setTitleX("#phi_t_r_a_j (deg)"); 
                    hi4.setTitleY("Counts"); 
                    hi4.setLineColor(k+1);
                    dgPhiT.addDataSet(hi4, i+type.getId()*3);
                }
                H2F hi5 = new H2F("hi"+type.getName()+region, "Region"+region, 80, -20, 20, 50, 0, 50);  
                hi5.setTitleX("edge (cm)"); 
                hi5.setTitleY("#theta_t_r_a_j (deg)"); 
                for(int k=0; k<6; k++) {
                    H1F hi7 = new H1F("hi"+type.getName()+region+k, "Region"+region, 100, -30, 30);  
                    hi7.setTitleX("#phi_t_r_a_j (deg)"); 
                    hi7.setTitleY("Counts"); 
                    hi7.setLineColor(k+1);
                    dgSecT.addDataSet(hi7, i+type.getId()*3);
                }
                dgTraj.addDataSet(hi1, i+type.getId()*3);
                dgPhi.addDataSet(hi2,  i+type.getId()*3);
                dgCut.addDataSet(hi3, i+type.getId()*3);
                dgEdge.addDataSet(hi5, i+type.getId()*3);
            }
        }
        if(!dgs.containsKey(TabGroup.KINEMATICS))
            dgs.put(TabGroup.KINEMATICS, new LinkedHashMap());
        dgs.get(TabGroup.KINEMATICS).put(1, dgEB);
        if(!dgs.containsKey(TabGroup.DC))
            dgs.put(TabGroup.DC, new LinkedHashMap());
        dgs.get(TabGroup.DC).put(1, dgTraj);
        dgs.get(TabGroup.DC).put(2, dgCut);
        dgs.get(TabGroup.DC).put(3, dgPhi);
        dgs.get(TabGroup.DC).put(4, dgEdge);
        dgs.get(TabGroup.DC).put(5, dgPhiT);
        dgs.get(TabGroup.DC).put(61, dgEdg1);
        dgs.get(TabGroup.DC).put(62, dgEdg2);
        dgs.get(TabGroup.DC).put(63, dgEdg3);
        dgs.get(TabGroup.DC).put(7, dgSecT);
        dgs.get(TabGroup.DC).put(8, dgEff);
    }
    
    
    public void processEvent(DataEvent event, int id, int particleId, int sec) {
        
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
//                if(pid==0) pid=charge*211;
                
                if(pindex==0) continue;
                if(pid==0 || clas!=2) continue;
                if(particleId!=0 && pid!=particleId) continue;
                
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

                if(sec!=0 && sector!=sec) continue;
                
                Particle part = new Particle(pid,
                                             recPart.getFloat("px", pindex),
                                             recPart.getFloat("py", pindex),
                                             recPart.getFloat("pz", pindex),
                                             recPart.getFloat("vx", pindex),
                                             recPart.getFloat("vy", pindex),
                                             recPart.getFloat("vz", pindex));
                if(part.p()<2) continue;
                
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
                    
                    if((region<3 && charge<0) || (region==3 && charge>0)) continue;
                    
                    //System.out.println(region);
                    dgs.get(TabGroup.DC).get(1).getH2F("hi"+type+region).fill(trajLocal.x(), trajLocal.y());
                    if(edge>delta) dgs.get(TabGroup.DC).get(2).getH2F("hi"+type+region).fill(trajLocal.x(), trajLocal.y());
                    dgs.get(TabGroup.DC).get(3).getH2F("hi"+type+region).fill(phi, theta);
                    if(edge>0) {
                        int side = (int) Math.signum(trajLocal.x());
                        dgs.get(TabGroup.DC).get(4).getH2F("hi"+type+region).fill(edge*side, theta);
                    }
                    if(Math.abs(theta-thetas[1])<dtheta){
                        for(int k=0; k<4; k++) {
                            if(k==0 || edge>-delta+delta*k) 
                                dgs.get(TabGroup.DC).get(5).getH1F("hi"+type+region+k).fill(phi);
                        }
                        if(edge>3) dgs.get(TabGroup.DC).get(7).getH1F("hi"+type+region+(sector-1)).fill(phi);
                    }
                }
                else if(detector==DetectorType.HTCC.getDetectorId()) {
                    dgs.get(TabGroup.KINEMATICS).get(1).getH2F("hiptheta"+type).fill(part.p(), Math.toDegrees(part.theta()));
                    dgs.get(TabGroup.KINEMATICS).get(1).getH2F("hiphitheta"+type).fill(Math.toDegrees(part.phi()), Math.toDegrees(part.theta()));
                }
            }		
        }

    }
    
    
    public void analyzeHisto(int nevents) {
        
        for(int id=0; id<2; id++) {
            String type = DataType.getType(id).getName();
            int icol = (id+1)*2;
            for(int region=1; region<4; region++) {
                H2F h2 = dgs.get(TabGroup.DC).get(4).getH2F("hi"+type+region);
                GraphErrors meanP = new GraphErrors("meanP");
                meanP.setTitleX("#theta (deg)");
                meanP.setTitleY("edge (cm)");
                meanP.setMarkerStyle(2);
                meanP.setMarkerSize(5);
                meanP.setMarkerColor(icol);
                GraphErrors meanN = new GraphErrors("meanN");
                meanN.setTitleX("#theta (deg)");
                meanN.setTitleY("edge (cm)");
                meanN.setMarkerStyle(4);
                meanN.setMarkerSize(5);
                meanN.setMarkerColor(icol);
                GraphErrors sigmaP = new GraphErrors();
                sigmaP.setTitleX("#theta (deg)");
                sigmaP.setTitleY("#sigma(edge) (cm)");
                sigmaP.setMarkerStyle(2);
                sigmaP.setMarkerSize(5);
                sigmaP.setMarkerColor(icol);
                GraphErrors sigmaN = new GraphErrors();
                sigmaN.setTitleX("#theta (deg)");
                sigmaN.setTitleY("#sigma(edge) (cm)");
                sigmaN.setMarkerStyle(4);
                sigmaN.setMarkerSize(5);
                sigmaN.setMarkerColor(icol);
                List<H1F> slices = h2.getSlicesY();
                for(int itheta=0; itheta<slices.size(); itheta++) {
                    H1F h1 = slices.get(itheta);
                    h1.setLineColor(icol);
                    h1.setTitle("#theta="+h2.getDataY(itheta)+"(deg)");
                    dgs.get(TabGroup.DC).get(60+region).addDataSet(h1, itheta);
                    if(Fiducial.getIntegralIDataSet(h1, 7.5, 8.5)<5) continue;
                    if(Fiducial.getIntegralIDataSet(h1, 9.5, 10.5) < Fiducial.getIntegralIDataSet(h1, 7.5, 8.5)*0.8) continue;
                    EdgeSlice f1 = new EdgeSlice("f1", -8, 8);
                    f1.setParameter(0, h1.getMax());
                    f1.setParameter(1, 3);
                    f1.setParameter(2, 1);
                    f1.setParameter(3, 0);
                    f1.setParameter(4, h1.getMax());
                    f1.setParameter(5,-3);
                    f1.setParameter(6, 1);
                    f1.setParameter(7, 0);
                    f1.setParLimits(2,   0.0,  5.0);
                    f1.setParLimits(3, -20.0, 20.0);
                    f1.setParLimits(6,   0.0,  5.0);
                    f1.setParLimits(7, -20.0, 20.0);
                    f1.setLineWidth(2);
                    DataFitter.fit(f1, h1, "Q");
                    if(f1.isFitValid() || true) {
                        meanP.addPoint(h2.getDataY(itheta), f1.getParameter(1), 0, f1.parameter(1).error());
                        meanN.addPoint(h2.getDataY(itheta),-f1.getParameter(5), 0, f1.parameter(5).error());
                        sigmaP.addPoint(h2.getDataY(itheta), Math.abs(f1.getParameter(2)), 0, f1.parameter(2).error());
                        sigmaN.addPoint(h2.getDataY(itheta), Math.abs(f1.getParameter(6)), 0, f1.parameter(6).error());
                    }
                }
                if(meanP.getDataSize(0)>0) {
                    dgs.get(TabGroup.DC).get(8).addDataSet(meanP, region-1);
                    dgs.get(TabGroup.DC).get(8).addDataSet(meanN, region-1);
                    dgs.get(TabGroup.DC).get(8).addDataSet(sigmaP, region-1 + 3);
                    dgs.get(TabGroup.DC).get(8).addDataSet(sigmaN, region-1 + 3);
                }
            }
        }
    }
    
    public class EdgeSlice extends Func1D{

            public EdgeSlice(String name, double min, double max) {
                    super(name, min, max);
                    super.addParameter("amp1");
                    super.addParameter("mean1");
                    super.addParameter("sigma1");
                    super.addParameter("slope1");
                    super.addParameter("amp2");
                    super.addParameter("mean2");
                    super.addParameter("sigma2");
                    super.addParameter("slope2");
            }

            @Override
            public double evaluate(double x){
                    double value = 0.0;
                    
                    if(x>0) {
                        if(x<this.getParameter(1)) 
                            value = this.getParameter(0)*FunctionFactory.gauss(x,this.getParameter(1), this.getParameter(2));
                        else
                            value = this.getParameter(0)+x*this.getParameter(3);
                    }
                    else {
                        if(x>this.getParameter(5)) 
                            value = this.getParameter(4)*FunctionFactory.gauss(x,this.getParameter(5), this.getParameter(6));
                        else
                            value = this.getParameter(4)+x*this.getParameter(7);
                    }
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
        
    
    public EmbeddedCanvasTabbed getCanvas() {
        EmbeddedCanvasTabbed canvas = null;
        for(TabGroup group : dgs.keySet()) {
            for(int tab : dgs.get(group).keySet()) {
                if(canvas == null)
                    canvas = new EmbeddedCanvasTabbed(group.getName()+tab);
                else
                    canvas.addCanvas(group.getName()+tab);
                canvas.getCanvas(group.getName()+tab).draw(dgs.get(group).get(tab));
                for(EmbeddedPad pad : canvas.getCanvas(group.getName()+tab).getCanvasPads()) {
                    pad.setTitleFont("Arial");
                    pad.setTitleFontSize(18);
                    if(!pad.getDatasetPlotters().isEmpty() && pad.getDatasetPlotters().get(0).getDataSet() instanceof H2F) 
                        pad.getAxisZ().setLog(true);
                    if(!pad.getDatasetPlotters().isEmpty() && pad.getDatasetPlotters().get(0).getDataSet() instanceof GraphErrors) {
                        pad.getAxisX().setRange(5, 36);
                        pad.getAxisY().setRange(0, 5);
                    }
                }
                if(group==TabGroup.DC && tab<3) this.drawDC(canvas.getCanvas(group.getName()+tab));
            }
        }
        return canvas;
    }
    
    private void drawDC(EmbeddedCanvas canvas) {
        
        for(int layer : new int[]{6, 18, 36}) {
            int superlayer = layer/6;
            for(int iedge=0; iedge<2; iedge++) {
                canvas.cd((superlayer-1)/2+iedge*3);
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
        parser.addOption("-pid","0", "particlr id (0=all)");
        parser.addOption("-s","0", "sector id (0=all)");
//        parser.addOption("-w", "1", "open graphical window (1) or run in batch mode (0)");
        parser.parse(args);
        
        
        int     maxEvents = parser.getOption("-n").intValue();
        int     pid       = parser.getOption("-pid").intValue();
        int     sector    = parser.getOption("-s").intValue();
        boolean window    = true;//parser.getOption("-w").intValue()==1;
        
        if(!window) System.setProperty("java.awt.headless", "true");
        DefaultLogger.debug();
        
        Fiducial analysis = new Fiducial();        

        List<String> inputFiles = parser.getInputList();
        
        ProgressPrintout progress = new ProgressPrintout();


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
        analysis.analyzeHisto(0);
    
        if(window) {
            JFrame frame = new JFrame("Fiducials");
            frame.setSize(1500,1000);
            frame.add(analysis.getCanvas());
            frame.setLocationRelativeTo(null);
            frame.setVisible(true);     
        }
    }
}



