#!/usr/bin/python


import sys
import getopt
import os
import xml.etree.ElementTree as ET
import logging

logging.basicConfig(filename='example.log',format = "%(levelname) -10s %(module)s:%(lineno)s %(message)s", level = logging.DEBUG)

class FlowItem:
    def __init__(self, stagename, input_entries, output_entries):
        self.stagename = stagename
        #These are lists of dicts. The object tag signifies the actual in the file system and is used for mapping 
        #the dataflow, whereas the port is used in binding to the tool
        self.input_entries = input_entries
        self.output_entries = output_entries
        
    def validate(self):
        #Checking if output names and input names do not occur more than once
        #TODO
        return True
        
    def portpresent(self, isInput, portname):
        if isInput:
            for input_entry in self.input_entries:
                if input_entry["port"] == portname:
                    return True
            return False
        else:
            for output_entry in self.output_entries:
                if output_entry["port"] == portname:
                    return True
            return False

class ToolItem:
    def __init__(self, toolname, toolpath, input_entries, output_entries):
        self.toolname = toolname
        self.input_entries = input_entries
        self.output_entries = output_entries
        self.toolpath = toolpath
        
    def validate(self):
        #TODO
        return True;
        
    def parameterpresent(self, isInput, parametername):
        
        if isInput:
            for input_entry in self.input_entries:
                if input_entry["name"] == parametername:
                    return True
            return False
        else:
            for output_entry in self.output_entries:
                if output_entry["name"] == parametername:
                    return True
            return False
            
            
class BindingItem:
    def __init__(self, flowname, toolname, input_entries, output_entries):
        self.toolname = toolname
        self.flowname = flowname
        self.input_entries = input_entries
        self.output_entries = output_entries
        
    def validate(self):
        #TODO
        return True;

class Workflow:
    def __init__(self, flow_xml_filename, binding_xml_name, tool_xml_name):
        self.error_list = []
        self.exemptflownames = ["begin", "end"]
        self.flows_list = self.parseflow(flow_xml_filename)
        self.tools_list = self.parsetool(tool_xml_name)
        self.binding_list = self.parseBinding(binding_xml_name)
        
        
        
        self.createaccessmaps()
        
    def printerrors(self):
        print "======================Error List=============================="
        logging.debug("======================Error List==============================")
        for error_item in self.error_list:
            print error_item
            logging.debug(error_item)
        
    #Creating easy lookups for binding objects
    def createaccessmaps(self):
        #Bindings
        binding_tool_map = {}
        binding_flow_map = {}
        for binding in self.binding_list:
            binding_tool_map[binding.toolname] = binding
            binding_flow_map[binding.flowname] = binding
        
        print binding_flow_map
        
        self.binding_tool_map = binding_tool_map
        self.binding_flow_map = binding_flow_map
        
        #Flows
        flow_map = {}
        for flow in self.flows_list:
            flow_map[flow.stagename] = flow
        self.flow_map = flow_map
        
        #tools
        tool_map = {}
        for tool in self.tools_list:
            tool_map[tool.toolname] = tool
        self.tool_map = tool_map
        
        
    #Parsing the Flow.xml
    def parseflow(self, flow_xml_filename):
        flows_list = []
        
        tree = ET.parse(flow_xml_filename)
        root = tree.getroot()
        for child in root:
            #Looking for only actions
            if(child.tag == "action"):
                print child.tag + "\t" + str(child.attrib)
                
                #Constructing parameters for flow
                stage_name = child.attrib["name"]
                input_entries = []
                output_entries = []
                
                #we can now create input and output for each flow object
                for dataflow in child:
                    print dataflow.tag
                    if dataflow.tag == "input":
                        input_entries.append(dataflow.attrib)
                        print dataflow.attrib
                    if dataflow.tag == "output":
                        output_entries.append(dataflow.attrib)
                        print dataflow.attrib
                
                flow_item = FlowItem(stage_name, input_entries, output_entries)
                flows_list.append(flow_item)
                
                
        return flows_list
        
    #Parsing the Tool.xml
    def parsetool(self, tool_xml_name):
        tools_list = []
        
        tree = ET.parse(tool_xml_name)
        root = tree.getroot()
        
        #First Getting All Tool Paths
        tool_path_present = {}
        for child in root:
            print child.tag + "\t" + str(child.attrib)
            if(child.tag == "pathSet"):
                for toolPathItem in child:
                    if toolPathItem.tag == "toolPath":
                        toolname = toolPathItem.attrib["tool"]
                        if toolname in tool_path_present:
                            print "Tool Redefinition";
                            self.error_list.append("Tool Path Redefinition in Tool.xml: " + toolname)
                        else:
                            tool_path_present[toolname] = child.attrib["base"] + "/" + toolPathItem.attrib["path"]
        
        print tool_path_present
        
        #Getting actual tool input and outputs
        for child in root:
            if(child.tag == "tool"):
                toolname = child.attrib["name"]
                tool_input = []
                tool_output = []
                
                for tool_parameters in child:
                    if tool_parameters.tag == "require":
                        tool_input.append(tool_parameters.attrib)
                    if tool_parameters.tag == "produce":
                        tool_output.append(tool_parameters.attrib)
                
                
                toolpath = ""
                #Checking if tool has path
                if toolname in tool_path_present:
                    toolpath = tool_path_present[toolname]
                else:
                    self.error_list.append("Missing path for tool: " + toolname)
                
                tool = ToolItem(toolname, toolpath, tool_input, tool_output)
                
                tools_list.append(tool)
                
        return tools_list
        
    #Parsing the binding xml file
    def parseBinding(self, binding_xml_name):
        binding_list = []
        
        tree = ET.parse(binding_xml_name)
        root = tree.getroot()
        
        for child in root :
            print child.tag + "\t" + str(child.attrib)
            if child.tag == "bind":
                flow_name = child.attrib["action"]
                
                if flow_name in self.exemptflownames:
                    continue
                
                tool_name = child.attrib["tool"]
                
                input_list = []
                output_list = []
                for binding_parameters in child:
                    if binding_parameters.tag == "inputAsRequirement":
                        input_list.append(binding_parameters.attrib)
                    if binding_parameters.tag == "productionToOutput":
                        output_list.append(binding_parameters.attrib)
                        
                bindingitem = BindingItem(flow_name, tool_name, input_list, output_list)
                binding_list.append(bindingitem)
                
        
        return binding_list
        
    @staticmethod
    def validate_flow_to_binding(flow_item, binding_item):
        output_errors = []
        
        #Checking the input entries
        binding_input_list = binding_item.input_entries
        binding_output_list = binding_item.output_entries
        
        #Checking Binding Import
        for binding_input in binding_input_list:
            binding_port_value = binding_input["port"]
            #Check if this is present
            if flow_item.portpresent(True, binding_port_value):
                print "PORT FOUND IN INPUT FLOW "  + flow_item.stagename
            else:
                output_errors.append("Port in binding not found in flow: " + flow_item.stagename)
        
        #Checking Binding Output
        print "CHECKING OUTPUT"
        
        for binding_output in binding_output_list:
            binding_port_value = binding_output["port"]
            #Check if this is present
            if flow_item.portpresent(False, binding_port_value):
                print "PORT FOUND IN OUTPUT FLOW " + flow_item.stagename
            else:
                output_errors.append("Port in binding not found in flow: " + flow_item.stagename)
        return output_errors
    
    @staticmethod
    def validate_binding_to_tool(binding_item, tool_item):
        output_errors = []
        
        #Checking from binding to tool
        binding_input_list = binding_item.input_entries
        binding_output_list = binding_item.output_entries
        
        
        for binding_input in binding_input_list:
            binding_requirement_name = binding_input["requirement"]
            #Check if this is present
            if tool_item.parameterpresent(True, binding_requirement_name):
                print "INPUT FOUND IN INPUT Tool "  + tool_item.toolname + " " + binding_requirement_name
            else:
                output_errors.append("Tool parameter in binding: " + binding_requirement_name + " not found in tool: " + tool_item.toolname)
        
        
        #Checking from tool to binding
        #TODO
        
        return output_errors
        
    #Validating the Tools, requires that the acccess maps are created
    def validate(self):
        #First validate Flow
        for flow_item in self.flows_list:
            if(flow_item.validate() == False):
                self.error_list.append("Flow Validation Error: " + flow_item.stagename)
        
        #Validating Tool
        for tool_item in self.tools_list:
            if(tool_item.validate() == False):
                self.error_list.append("Tool Validation Error: " + tool_item.toolname)
                
        #Validating no cycles in data flow
        #We'll punt on this for a while, we'll only check that it doesn't go from input to output on the same node

                
        #Validating the connection from tool to flow through binding
        for flow in self.flows_list:
            flowname = flow.stagename
            
            if flowname in self.exemptflownames:
                continue
            
            if flowname in self.binding_flow_map:
                binding_item = self.binding_flow_map[flowname]
                self.error_list.extend(Workflow.validate_flow_to_binding(flow, binding_item))
                
                #Validate Binding to Tool
                tool_name = binding_item.toolname
                
                if tool_name in self.tool_map:
                    self.error_list.extend(Workflow.validate_binding_to_tool(binding_item, self.tool_map[tool_name]))
                else:
                    self.error_list.append("Tool not found for flow: " + flowname + " should be named: " + tool_name)
                
            else:
                self.error_list.append("Binding not found for flow: " + flowname)
            
            
def usage():
    print "<input flow.xml> <input binding.xml> <input tool.xml>"
    


    
def main():
    usage()
    
    workflow = Workflow(sys.argv[1], sys.argv[2], sys.argv[3])
    workflow.validate() == True
    
    workflow.printerrors();
    
if __name__ == "__main__":
    main()
