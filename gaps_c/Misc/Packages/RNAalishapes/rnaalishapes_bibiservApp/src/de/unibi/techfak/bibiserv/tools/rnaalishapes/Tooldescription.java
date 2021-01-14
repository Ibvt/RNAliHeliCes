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
 * Contributor(s): Jan Krueger - jkrueger(at)cebitec.uni-bielefeld.de
 * 
 */
package 

de.unibi.techfak.bibiserv.tools.rnaalishapes;

import de.unibi.techfak.bibiserv.cms.Tfunction;
import de.unibi.techfak.bibiserv.cms.TrunnableItem;
import java.util.HashMap;
import java.util.Map;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBElement;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Unmarshaller;
import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import org.apache.log4j.Logger;
import org.springframework.beans.factory.InitializingBean;
import org.w3c.dom.Element;

/**
 * Bean class which returns tools current tool description
 *
 *
 * @author Jan Krueger - jkrueger(at)cebitec.uni-bielefeld.de
 */
public class Tooldescription implements InitializingBean {

    private static Logger log = Logger.getLogger(Tooldescription.class);
    private TrunnableItem tooldescription;
    private Map<String, Tfunction> functionmap;
    private Element e_tooldescription;

    /**
     * Service method that return a tool description. The tool description is
     * currently an RunableItem Object.
     *
     * @return Return the runnableItem.
     */
    public TrunnableItem getToolDescription() {
        return tooldescription;
    }

    /**
     * Return a function object belongig to given id.
     *
     * @param id
     * @return
     */
    public Tfunction getFunction(String id) {
        return functionmap.get(id);
    }

    /**
     * Service method that returns a tool description as DOM element.
     *
     *
     * @return
     */
    public Element getToolDescriptionAsDOM() {
        if (e_tooldescription == null) {
            try {
                /*
                 * read param from (xml-file)
                 */
                DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
                dbf.setNamespaceAware(true);
                DocumentBuilder db = dbf.newDocumentBuilder();
                /*
                 * read description from file
                 */
                e_tooldescription = db.parse(getClass().getResourceAsStream("/runnableitem.xml")).getDocumentElement();
            } catch (Exception e) {
                throw new RuntimeException(e);
            }
        }
        return e_tooldescription;
    }

    @Override
    public void afterPropertiesSet() throws Exception {
        if (tooldescription == null) {
            try {
                // create JAXB Context
                JAXBContext jaxbcontext = JAXBContext.newInstance(TrunnableItem.class);
                // create Unmarshaeller
                Unmarshaller unmarshaller = jaxbcontext.createUnmarshaller();
                // load runableitem from resource and ...
                JAXBElement<TrunnableItem> jaxbe = (JAXBElement) unmarshaller.unmarshal(getClass().getResourceAsStream("/runnableitem.xml"));
                // ... set tooldescription
                tooldescription = jaxbe.getValue();
                // store all functions in function map for faster access
                functionmap = new HashMap<String, Tfunction>();
                if (tooldescription.isSetExecutable()) {
                    for (Tfunction tf : tooldescription.getExecutable().getFunction()) {
                        functionmap.put(tf.getId(), tf);
                    }
                }
            } catch (JAXBException e) {
                log.fatal("A JAXBException occurred while reading RunnableItem as Stream from Resource! Errormessage was : " + e.getMessage());
                throw new RuntimeException("A JAXBException occurred while reading RunnableItem as Stream from Resource! Errormessage was : " + e.getMessage(), e);
            }
        }
    }
}
