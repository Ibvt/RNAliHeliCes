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
 * "Portions Copyrighted [year] [name of copyright owner]"
 *
 * Contributor(s):
 *
 */
/**
 * Implementation class  for function rnaalishapes_function_mfe_threadworker, implements all request and
 * response methods ...
 *
 * <b>Attention: This is autogenerated code. </b>
 *
 *
 * @author Thomas Gatter - tgatter[aet]cebitec.uni-bielefeld.de (template)
 */
package de.unibi.techfak.bibiserv.tools.rnaalishapes;

import de.unibi.techfak.bibiserv.util.ontoaccess.bibiontotypes.OntoRepresentation;
import de.unibi.techfak.bibiserv.util.ontoaccess.bibiontotypes.OntoRepresentation.datastructure;
import de.unibi.cebitec.bibiserv.utils.connect.AWSValidationConnection;
import de.unibi.cebitec.bibiserv.utils.connect.URLValidationConnection;
import de.unibi.techfak.bibiserv.BiBiTools;
import de.unibi.techfak.bibiserv.CmdLineInfo;
import de.unibi.cebitec.bibiserv.utils.UniversalValidator;
import de.unibi.cebitec.bibiserv.utils.UniversalConverter;
import de.unibi.techfak.bibiserv.Call;
import de.unibi.techfak.bibiserv.CallFactory;
import de.unibi.techfak.bibiserv.Status;
import de.unibi.techfak.bibiserv.exception.BiBiToolsException;
import de.unibi.techfak.bibiserv.exception.DBConnectionException;
import de.unibi.techfak.bibiserv.exception.IdNotFoundException;
import de.unibi.techfak.bibiserv.util.LogShed;
import de.unibi.techfak.bibiserv.util.Pair;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.List;

import org.apache.log4j.Logger;

public class rnaalishapes_function_mfe_threadworker {


private String fct_id = "rnaalishapes_function_mfe";


private static Logger log = Logger.getLogger(rnaalishapes_function_mfe_threadworker.class);

public void startRequestThread(final BiBiTools bibitools, final List<Pair<String,String>> paramlist, final String accesskey, final String secretkey, final String uploadbucket, final String uploadfolder, final File out,final Object input_0,final OntoRepresentation representation_0,final OntoRepresentation representation_input_0,final boolean streamsSupported_0,final boolean skipValidation_0) throws BiBiToolsException {


        final Status status = bibitools.getStatus();

        // create new ProcessingThread ...
        Thread processing = new Thread(new Runnable() {

            @Override
            public void run() {
                try {
                    try {
                        // parse and check parameter
                        HashMap<String,String> paramhash = bibitools.checkAndParseParam(paramlist, fct_id, new LogShed());
                        // check and parse input
                        HashMap<String,String> inputhash = new HashMap<String,String>();
                        String prefix = "";
                        String postfix = "";
                        CmdLineInfo generatedInfo = new CmdLineInfo();
                        
if (input_0 instanceof AWSValidationConnection) {
	String validatorImpl_0;
	List<String> converterChain_0;
	if(skipValidation_0) {
		validatorImpl_0 = null;
		converterChain_0 = null;
	} else {
		validatorImpl_0 = UniversalValidator.getValidatorImlementation(representation_0);
		converterChain_0 = UniversalConverter.getStreamConverterOrder(representation_0, representation_input_0);
	}

	AWSValidationConnection connection = (AWSValidationConnection) input_0;

	postfix += bibitools.parseInputAWS("rnaalishapes_input_rna_sequence_alignment",inputhash, connection.getBucket(), connection.getFile(), accesskey, secretkey, generatedInfo, validatorImpl_0, converterChain_0, (representation_0.getContent()!=null) ? representation_0.getContent().name() : "NULL" , (representation_0.getStrictness()!=null) ? representation_0.getStrictness().name() : "NULL" , (representation_0.getCardinality()!=null) ? representation_0.getCardinality().name() : "NULL" , representation_0.getStructure().equals(datastructure.ALIGNMENT) || representation_0.getStructure().equals(datastructure.STRUCTUREALIGNMENT), streamsSupported_0);
} else if (input_0 instanceof URLValidationConnection) {
	String validatorImpl_0;
	List<String> converterChain_0;
	if(skipValidation_0) {
		validatorImpl_0 = null;
		converterChain_0 = null;
	} else {
		validatorImpl_0 = UniversalValidator.getValidatorImlementation(representation_0);
		converterChain_0 = UniversalConverter.getStreamConverterOrder(representation_0, representation_input_0);
	}

	URLValidationConnection connection = (URLValidationConnection) input_0;
	postfix += bibitools.parseInputURL("rnaalishapes_input_rna_sequence_alignment",inputhash, connection.getUrl(), generatedInfo, validatorImpl_0, converterChain_0, (representation_0.getContent()!=null) ? representation_0.getContent().name() : "NULL" , (representation_0.getStrictness()!=null) ? representation_0.getStrictness().name() : "NULL" , (representation_0.getCardinality()!=null) ? representation_0.getCardinality().name() : "NULL" , representation_0.getStructure().equals(datastructure.ALIGNMENT) || representation_0.getStructure().equals(datastructure.STRUCTUREALIGNMENT), streamsSupported_0);
} else {
	postfix += bibitools.parseInput("rnaalishapes_input_rna_sequence_alignment",inputhash, input_0, representation_0.getType().name(), representation_0.getImplementationType());
}

                        // merge both hashmaps
                        inputhash.putAll(paramhash);

                        generateUploads(bibitools, status, out, uploadbucket, uploadfolder, accesskey, secretkey, generatedInfo);
                        
                        /*
                         * TG 01/13: Change behaviour of getOutputFile if needed!
                         */
                        String toolCmd = bibitools.generateCmdLineString(fct_id, inputhash, prefix, postfix);
                        String cmdline = bibitools.generateStreamCmdScriptString(toolCmd, generatedInfo, out);
                        
                        
                        // do some preprocessing
                        status.setStatuscode(602);
                        // do processing
                        CallFactory cf = CallFactory.newInstance();
                        Call call = cf.newCall(bibitools);
                        if (!call.call(cmdline)) {
                            // tool execution failed
                            String message = BiBiTools.i2s(new InputStreamReader(call.getStdErrStream()));
                            if(status.getStatuscode()<700) { // only set status if not already set
                                status.setStatuscode(703, message);
                            }
                            log.info("Tool exec failed: "+status.getStatuscode()+" "+status.getDescription()+". ErrorStream: "+message);
                            return;
                        }
                        // do some postprocessing
                        status.setStatuscode(605);
                        // processing finished
                        status.setStatuscode(600);
                    } catch (BiBiToolsException e) {
                        // setting of status is not necessary, because
                        // this should be done by the Call implementation
                         log.error("BiBiToolsException "+e);
                    } catch (IOException e) {
                        log.error("IOException "+e);
                        status.setStatuscode(720);
                    } catch (IdNotFoundException e){
                         status.setStatuscode(722, "Internal Resource Error");
                    }
                } catch (DBConnectionException e) {
                    log.fatal(e.getMessage(),e);
                }
            }
        });
        // ... and start it 
        processing.start();
    

}

	
   public void generateUploads(BiBiTools bibitools, Status status, File out, String uploadbucket, String uploadfolder, String accesskey, String secretkey, CmdLineInfo generatedInfo) throws BiBiToolsException, DBConnectionException {
        if (uploadbucket == null || uploadfolder == null) { // null: user did not specify upload!
            return; 
        }
        String folder;
        try {
            folder = bibitools.getSpoolDir().toURI().relativize(out.getParentFile().toURI()).getPath();
        } catch (FileNotFoundException e) {
            status.setStatuscode(701, "Could not find upload folder..");
            throw new BiBiToolsException(status.toString(), e);
        }
        bibitools.parseUpload(uploadbucket, uploadfolder, accesskey, secretkey, folder, out.getName(), generatedInfo);
        
        // additional output files
        
    }
}
