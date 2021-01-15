/*
 * DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS HEADER.
 *
 * Copyright 2010 BiBiServ Curator Team, http://bibiserv.cebitec.uni-bielefeld.de,
 * All rights reserved.
 *
 * The contents of this file are subject to the terms of the Common
 * Development and Distribution License("CDDL") (the "License"). You
 * may not use this file except in compliance with the License. You can
 * obtain a copy of the License at http://www.sun.com/cddl/cddl.html
 *
 * See the License for the specific language governing permissions and
 * limitations under the License.  When distributing the software, include
 * this License Header Notice in each file.  If applicable, add the following
 * below the License Header, with the fields enclosed by brackets [] replaced
 *  by your own identifying information:
 *
 * "Portions Copyrighted 2010 BiBiServ Curator Team, http://bibiserv.cebitec.uni-bielefeld.de"
 *
 * Contributor(s):  Daniel Hagemeier - dhagemei[aet]cebitec.uni-bielefeld.de
 *                  Armin Toepfer - atoepfer[aet]cebitec.uni-bielefeld.de
 *                  Jan Krueger - jkrueger[aet]cebitec.uni-bielefeld.de
 *
 */
package de.unibi.techfak.bibiserv.tools.rnaalishapes.web;

import de.unibi.cebitec.bibiserv.util.convert.ConversionException;
import de.unibi.cebitec.bibiserv.util.validate.ValidationException;
import de.unibi.techfak.bibiserv.util.ontoaccess.bibiontotypes.OntoRepresentation;
import de.unibi.cebitec.bibiserv.utils.SiblingGetter;
import de.unibi.cebitec.bibiserv.statistics.logging.StatsLoggerI;
import de.unibi.techfak.bibiserv.cms.Texample;
import de.unibi.techfak.bibiserv.cms.Texample.Prop;

import de.unibi.techfak.bibiserv.web.beans.session.MessagesInterface;
import de.unibi.techfak.bibiserv.util.Pair;
import de.unibi.techfak.bibiserv.exception.BiBiToolsException;
import de.unibi.techfak.bibiserv.tools.rnaalishapes.ParameterDependencies;
import de.unibi.techfak.bibiserv.tools.rnaalishapes.Tooldescription;
import de.unibi.techfak.bibiserv.tools.rnaalishapes.rnaalishapes_function_convert;
import de.unibi.techfak.bibiserv.util.dependencyparser.ConstraintHashMap;
import de.unibi.techfak.bibiserv.util.dependencyparser.Constraints;
import de.unibi.techfak.bibiserv.util.dependencyparser.DependencyException;
import de.unibi.techfak.bibiserv.util.dependencyparser.Id;
import de.unibi.techfak.bibiserv.util.dependencyparser.Node;
import de.unibi.techfak.bibiserv.util.dependencyparser.ParameterWrapper;
import de.unibi.techfak.bibiserv.web.beans.session.AbstractCloudInputBean;
import java.io.IOException;
import java.util.ArrayList;

import java.util.List;
import org.primefaces.context.RequestContext;
import javax.faces.context.FacesContext;
import javax.faces.event.ActionEvent;
import javax.faces.event.ValueChangeEvent;

import org.apache.log4j.Logger;
import org.springframework.beans.factory.DisposableBean;
import org.springframework.beans.factory.InitializingBean;

/**
 * This is a autogenerated bean controller template class for function <i>rnaalishapes_function_convert_controller</i>
 *
 * @author Daniel Hagemeier - dhagemei[aet]cebitec.uni-bielefeld.de
 *         Armin Toepfer - atoepfer[aet]cebitec.uni-bielefeld.de
 *         Jan Krueger - jkrueger[aet]cebitec.uni-bielefeld.de
 */
public class rnaalishapes_function_convert_controller implements InitializingBean, DisposableBean {

    private final static Logger log = Logger.getLogger(de.unibi.techfak.bibiserv.tools.rnaalishapes.web.rnaalishapes_function_convert_controller.class);
    private final String br = System.getProperty("line.separator");



        /* #########################################
     * #   DI method setguugle and field   #
     * ######################################### */
    public rnaalishapes_function function;

    public void setFunction(rnaalishapes_function function) {
        this.function = function;
    }

    /* #########################################
     * #     DI method setMessages             #
     * ######################################### */
    private MessagesInterface messages;

    public void setMessages(MessagesInterface messages) {
        this.messages = messages;
    }

    /* #########################################
     * # DI method setParameterDependencies    #
     * ######################################### */
    private ParameterDependencies pdp;

    public void setParameterDependencies(ParameterDependencies pdp) {
        this.pdp = pdp;
    }

    /* #########################################
     * #     DI method setStatsLogger          #
     * ######################################### */
    private StatsLoggerI statsLogger;

    public void setStatsLogger(StatsLoggerI statsLogger) {
        this.statsLogger = statsLogger;
    }

    
    
    /* #########################################
     * #     DI method setInput(s)             #
     * ######################################### */
    
    public AbstractCloudInputBean rnaalishapes_input_rna_secondary_structure;

    public void setInput(AbstractCloudInputBean rnaalishapes_input_rna_secondary_structure) {
        this.rnaalishapes_input_rna_secondary_structure = rnaalishapes_input_rna_secondary_structure;
    }
    

    /* #########################################
     * #     DI method setParam                 #
     * ######################################### */
    public rnaalishapes_function_convert_param param;

    public void setParam(rnaalishapes_function_convert_param param) {
        this.param = param;
    }

    /* #########################################
     * #     DI method setResult               #
     * ######################################### */
    public rnaalishapes_function_convert_result result;

    public void setResult(rnaalishapes_function_convert_result result) {
        this.result = result;
    }
    
    /* #########################################
     * #     DI method setResultHandler        #
     * ######################################### */
    public rnaalishapes_function_convert_resulthandler resulthandler;

    public void setResulthandler(rnaalishapes_function_convert_resulthandler resulthandler) {
        this.resulthandler = resulthandler;
    }

    /* #########################################
     * #     DI method setExecfunction          #
     * ######################################### */
    public rnaalishapes_function_convert execfunction;

    public void setExecfunction(rnaalishapes_function_convert execfunction) {
        this.execfunction = execfunction;
    }
    
      /* #########################################
     * #     DI method setTooldescription       #
     * ######################################### */
    public Tooldescription tooldescription;

    public void setTooldescription(Tooldescription tooldescription) {
        this.tooldescription = tooldescription;
    }
    
    /* #########################################
     * # Returns all Representations as String #
     * ######################################### */
    
    public String getSupportedFormats() {
        OntoRepresentation base = execfunction.getRepresentationInput();
        
         List<OntoRepresentation> possibleRepresentations = new ArrayList<OntoRepresentation>();
        possibleRepresentations.add(base);
        possibleRepresentations.addAll(SiblingGetter.getSiblingsConvertableTo(base));

        StringBuilder sup = new StringBuilder();
        boolean first = true;
        for (OntoRepresentation representation : possibleRepresentations) {
            if(first){
                 sup.append(representation.getLabel());
                 first = false;
             } else {
                 sup.append(", ").append(representation.getLabel());
             }
        }
        return sup.toString();
    }
    
    public String getSupportedFormatsStream() {
        OntoRepresentation base = execfunction.getRepresentationInput();
        
         List<OntoRepresentation> possibleRepresentations = new ArrayList<OntoRepresentation>();
        possibleRepresentations.add(base);
        possibleRepresentations.addAll(SiblingGetter.getSiblingsStreamConvertableTo(base));

        StringBuilder sup = new StringBuilder();
        boolean first = true;
        for (OntoRepresentation representation : possibleRepresentations) {
            if(first){
                 sup.append(representation.getLabel());
                 first = false;
             } else {
                 sup.append(", ").append(representation.getLabel());
             }
        }
        return sup.toString();
    }
    
    
    public String getToolMainFormat() {
        OntoRepresentation base = execfunction.getRepresentationInput();
        return base.getLabel();
    }
    
    
    
    /* ################################################
     * # Implementation of Interface InitializingBean #
     * ################################################ */
   
    private Node dptree;
    private ParameterWrapper parameterwrapper;

    public void afterPropertiesSet() throws Exception {
        // create new parameterwrapper
        parameterwrapper = new ParameterWrapper();
        // get DependencyTree
        dptree = pdp.getDependencyTree("rnaalishapes_function_convert", parameterwrapper);
    }

    /* ###############################################
     * # Implementation of Interface Disposable Bean #
     * ############################################### */
    public void destroy() throws Exception {       
    }

    /* #########################################
     * #     waiting dialig showing beans      #
     * ######################################### */
    

    public boolean isInputValidated() {
        return true 
                
                && rnaalishapes_input_rna_secondary_structure.isValid()
                ;
                
    }
    
   
    
    
   /* ###############################################################################
     * # Callbackfunction sending success of validation to client with RequestContext#
     * ############################################################################## */
     public void callback() {
      
        RequestContext context = RequestContext.getCurrentInstance();
        context.addCallbackParam("validated",true&&rnaalishapes_input_rna_secondary_structure.isValidated());
                 
         if (true&&rnaalishapes_input_rna_secondary_structure.isValid()) {     
             log.debug("Request to callback validation....");
             context.addCallbackParam("redirection",true ); 
             context.addCallbackParam("location", "2"); 
             // resetting values for the user returning to input page...
             
             
             rnaalishapes_input_rna_secondary_structure.setValidated(false);
             if(rnaalishapes_input_rna_secondary_structure.getInput().getRepresentations().size()>1) {
                 context.addCallbackParam("location", "1b"); 
             }
             
            }

    }
    
   
    /* #########################################
     * # Validation method for this function   #
     * # initialised by Calculate-Button       #
     * ######################################### */
    public void validate(ActionEvent e) {
        log.info("call validate ...");
        

        // set needed i18n properties in string array, facescontext.currentInstance return
        // null value, when called from new internal thread
        final String[] args = new String[]{
            messages.property("de.unibi.techfak.bibiserv.bibimainapp.input.UPLOAD"),
            messages.property("de.unibi.techfak.bibiserv.bibimainapp.input.NOFILE"),
            messages.property("de.unibi.techfak.bibiserv.bibimainapp.input.COPYPASTE"),
            messages.property("de.unibi.techfak.bibiserv.bibimainapp.input.NODATA"),
            messages.property("de.unibi.techfak.bibiserv.bibimainapp.input.NOXML"),
            messages.property("de.unibi.techfak.bibiserv.bibimainapp.input.CHECK"),
            messages.property("de.unibi.techfak.bibiserv.bibimainapp.input.MESSAGE")
        };

        // validation can be time consuming, so start validation in an extra thread
        Thread validation = new Thread(new Runnable() {

            public void run() {
                OntoRepresentation target = null;
                
                target = execfunction.getRepresentationInput();
                rnaalishapes_input_rna_secondary_structure.determineSourceAndValidate(args,target);
                
                log.info("Validation finished ...");
               
                    //disable error field(s)
                    
                    if (rnaalishapes_input_rna_secondary_structure.isValid()){
                        rnaalishapes_input_rna_secondary_structure.setShowInfo(false);
                    }
                    
                    
                    
                } 
               
            
        });
        
        validation.start();
        
    }


    /* #########################################
     * # Calculation method for this function  #
     * # initialised by Calculate-Button       #
     * ######################################### */
    public void calculate(ActionEvent e) {
        
        // only resume if the resulthandler is valid
        if(!resulthandler.validate()){
            return;
        }
        
        try {
            log.info("Called calculate");

            final List<Pair<String, String>> pairlist = param.getParameterList();
            
            /* ---------------------- Stats ---------------------- */
            statsLogger.logSubmit("rnaalishapes_submission", function.getSessionId(), defaultParams(), defaultInputs());
            long startTime = System.currentTimeMillis();

            try {
                /* check each combination of inputs and call appropriate function */
                String requestid;
                String secretkey = resulthandler.getSecretKey();
                String accesskey = resulthandler.getAccessKey();
                if(resulthandler.getResultHandling()==resulthandler.getResultHandling().s3upload) {
                    String uploadbucket = resulthandler.getSelected_item_bucket();
                    String uploadfolder = resulthandler.getSubfolder();
                    result.setShowMessage(true);
                    result.setCalculateMessage(messages.property("de.unibi.techfak.bibiserv.bibimainapp.result.UPLOAD_COMPLETED"));
                    
                    requestid = execfunction.request(pairlist, accesskey, secretkey, uploadbucket, uploadfolder, rnaalishapes_input_rna_secondary_structure.getInput().getInput(), rnaalishapes_input_rna_secondary_structure.getInput().getChosen(), rnaalishapes_input_rna_secondary_structure.supportsStreamedInput(), rnaalishapes_input_rna_secondary_structure.getInput().isSkipValidation());

                } else {
                    result.setShowMessage(false);
                    requestid = execfunction.request(pairlist, accesskey, secretkey , rnaalishapes_input_rna_secondary_structure.getInput().getInput(), rnaalishapes_input_rna_secondary_structure.getInput().getChosen(), rnaalishapes_input_rna_secondary_structure.supportsStreamedInput(), rnaalishapes_input_rna_secondary_structure.getInput().isSkipValidation());

                }
                result.setBibiservid(requestid);

            } catch (BiBiToolsException ex) {
                log.fatal("BiBiToolException occurred! ", ex);
            } catch (ConversionException ex) {
                log.fatal("ConversionException occurred (should never happend at this point)! ", ex);
            } catch (ValidationException ex) {
                log.fatal("ValidateException occurred (should never happend at this point)! ", ex);
            }
            statsLogger.logRuntime("rnaalishapes_submission", function.getSessionId(), result.getBibiservid(), (System.currentTimeMillis() - startTime) / 1000, result.getStatuscode(), result.getStatusdescription());

            FacesContext.getCurrentInstance().getExternalContext().redirect("/rnaalishapes?viewType=submission&subType=rnaalishapes_function_convert_result");

        } catch (IOException ex){
            log.fatal(ex.getMessage());
        } 
        
    }
    
    
    /* #########################################
     * # Calculation method for this function  #
     * # initialised by Calculate-Button       #
     * ######################################### */
    public void paramSet(ActionEvent e) {
        log.info("Called Paramset");

        // if parameters aren't valid ,return without calculation
        if (!param.isValid()) {
            log.debug("Parameters aren't valid, abbort further calculation ...");
            return;
        }

        /* -------------------- Parameter -------------------- */
        final List<Pair<String, String>> pairlist = param.getParameterList();

        parameterwrapper.setParameter(pairlist);


        // if tree evaluates to true we can start the calculation
        try {
            if (dptree.evaluate()) {

                FacesContext.getCurrentInstance().getExternalContext().redirect("/rnaalishapes?viewType=submission&subType=rnaalishapes_function_convert_p_3");

            } else {
                ConstraintHashMap chm = dptree.getMissingConstraints();
                FacesContext fc = FacesContext.getCurrentInstance();
                for (Id id : chm.keySet()) {
                    List<Constraints> constraints = chm.get(id);
                    // Constraints list is null  or is empty parameter 'id' is missing without further constraints.
                    if (constraints == null || constraints.isEmpty()) {
                        param.addFaultmsg(id.getId()+"_rnaalishapes_function_convert","Parameter '"+id.getId()+"' is mandantory !");
                    } else {
                        StringBuilder sb = new StringBuilder("Parameter  '"+id.getId()+" is mandantory and must fullfill the constraints : ");
                        for (Constraints constraint : constraints) {
                            sb.append("'");
                            sb.append(constraint.toString());
                            sb.append("' ");
                        }
                        param.addFaultmsg(id.getId()+"_rnaalishapes_function_convert",sb.toString());
                    }
                }
            }
        } catch (DependencyException ex) {
            log.fatal(ex.getMessage(), ex);
        } catch (IOException ex){
            log.fatal(ex.getMessage());
        } 
    }

    /* ##################
     * # reset function #
     * ################## */
    public void reset() {
        resetParam(null);
        resetInput(null);
    }

    public void resetParam(ActionEvent e) {
        param.reset();
    }

    public void resetInput(ActionEvent e) {
        
        rnaalishapes_input_rna_secondary_structure.reset();
        
        
    }
    
      /* ######################
     * # example function(s) #
     * ####################### */
    
    

    /**  
     * Return the count of examples for this function 
     * 
     * @return number of examples for current function
     */
    public int getExampleCount(){             
        return tooldescription.getFunction("rnaalishapes_function_convert").getExample().size();
    }
       
    /**
     * ExampleActionListener, used by example button and example som chooser.
     * 
     * Replace values with the given example values.
     */
    public void example(int i) {
       Texample te = tooldescription.getFunction("rnaalishapes_function_convert").getExample().get(i);
       // iterate over all properties
       for (Prop p : te.getProp()) {
           rnaalishapes_input_rna_secondary_structure.checkAndSet(p);
           
           
           param.checkAndSet(p);
       }
    }
    
     public String exampleName(int i) {
        return tooldescription.getFunction("rnaalishapes_function_convert").getExample().get(i).getName().get(0).getValue();
    }
     
     /**
     * Return true in the case that input(s) follows (one of) the input example
     * 
     * @return 
     */
    
     public boolean defaultInputs(){
         return false;
     }
     
     /**
      * 
      * @return  Return true in the case the param(s) follows the example
      */
     public boolean defaultParams(){
         return false;
     }
}