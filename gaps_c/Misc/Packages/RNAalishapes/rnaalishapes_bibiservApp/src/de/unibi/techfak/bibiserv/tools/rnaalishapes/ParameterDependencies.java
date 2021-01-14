/*
 * DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS HEADER.
 * 
 * Copyright 2011 BiBiServ Curator Team, http://bibiserv.cebitec.uni-bielefeld.de, 
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

package de.unibi.techfak.bibiserv.tools.rnaalishapes;

import de.unibi.techfak.bibiserv.util.dependencyparser.DependencyParser;
import de.unibi.techfak.bibiserv.util.dependencyparser.Node;
import de.unibi.techfak.bibiserv.util.dependencyparser.ParameterWrapper;
import org.apache.log4j.Logger;
import org.springframework.beans.factory.InitializingBean;

/**
 * BeanWrapper for DependencyParser.
 *
 *
 * @author Jan Krueger - jkrueger(at)cebitec.uni-bielefeld.de
 */
public class ParameterDependencies implements  InitializingBean{

    private final static Logger log = Logger.getLogger(ParameterDependencies.class);

    private DependencyParser dp;

    private Tooldescription tooldesription;

    public void setTooldesription(Tooldescription tooldesription) {
        this.tooldesription = tooldesription;
    }

    public void afterPropertiesSet() throws Exception {
       dp = new DependencyParser();
       dp.setTooldescription(tooldesription.getToolDescriptionAsDOM());
    }

    public Node getDependencyTree(String function_id,ParameterWrapper pw) {
        try {
            dp.setFunctionId(function_id);
            dp.setParameterWrapper(pw);
            return dp.generate();
        } catch (Exception e){
            log.fatal("Exception occurred while generate DependencyTree.",e);
            throw new RuntimeException(e);
        }
    }

}
