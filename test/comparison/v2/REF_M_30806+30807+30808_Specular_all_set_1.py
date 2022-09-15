# Cross-section: Off_Off
# Run:30806
######################################################################
# Python Script Generated by GeneratePythonScript Algorithm
######################################################################
LoadEventNexus(Filename="/SNS/REF_M/IPTS-21391/nexus/REF_M_30794.nxs.h5", OutputWorkspace="raw_events")
FilterByLogValue(
    InputWorkspace="raw_events",
    OutputWorkspace="30794_entry-Off_Off",
    LogName="BL4A:SF:ICP:getDI",
    MinimumValue=15,
    MaximumValue=15,
    TimeTolerance=0.10000000000000001,
    LogBoundary="Left",
)
GroupWorkspaces(InputWorkspaces="30794_entry-Off_Off", OutputWorkspace="wsg")
LoadEventNexus(Filename="/SNS/REF_M/IPTS-21391/nexus/REF_M_30806.nxs.h5", OutputWorkspace="raw_events")
FilterByLogValue(
    InputWorkspace="raw_events",
    OutputWorkspace="30806_entry-Off_Off",
    LogName="BL4A:SF:ICP:getDI",
    MinimumValue=15,
    MaximumValue=15,
    TimeTolerance=0.10000000000000001,
    LogBoundary="Left",
)
GroupWorkspaces(InputWorkspaces="30806_entry-Off_Off", OutputWorkspace="wsg")
MagnetismReflectometryReduction(
    InputWorkspace="wsg",
    NormalizationWorkspace="30794_entry-Off_Off",
    SignalPeakPixelRange="181,195",
    SignalBackgroundPixelRange="49,88",
    NormPeakPixelRange="202,216",
    NormBackgroundPixelRange="94,104",
    LowResDataAxisPixelRange="69,172",
    LowResNormAxisPixelRange="71,175",
    TimeAxisRange="11413.6,45388.8",
    RoundUpPixel=False,
    SpecularPixel=187.30000000000001,
    FinalRebin=False,
    QMin=0.001,
    QStep=-0.02,
    TimeAxisStep=400,
    ConstQTrim=0.10000000000000001,
    OutputWorkspace="r30806",
)
Scale(InputWorkspace="r30806", OutputWorkspace="r30806", Factor=1.6068516739750773)
Scale(InputWorkspace="r30806", OutputWorkspace="r30806_scaled", Factor=2.4266100000000002)
AddSampleLog(Workspace="r30806_scaled", LogName="scaling_factor", LogText="2.42661", LogType="Number")
# Run:30807
######################################################################
# Python Script Generated by GeneratePythonScript Algorithm
######################################################################
LoadEventNexus(Filename="/SNS/REF_M/IPTS-21391/nexus/REF_M_30794.nxs.h5", OutputWorkspace="raw_events")
FilterByLogValue(
    InputWorkspace="raw_events",
    OutputWorkspace="30794_entry-Off_Off",
    LogName="BL4A:SF:ICP:getDI",
    MinimumValue=15,
    MaximumValue=15,
    TimeTolerance=0.10000000000000001,
    LogBoundary="Left",
)
GroupWorkspaces(InputWorkspaces="30794_entry-Off_Off", OutputWorkspace="wsg")
LoadEventNexus(Filename="/SNS/REF_M/IPTS-21391/nexus/REF_M_30807.nxs.h5", OutputWorkspace="raw_events")
FilterByLogValue(
    InputWorkspace="raw_events",
    OutputWorkspace="30807_entry-Off_Off",
    LogName="BL4A:SF:ICP:getDI",
    MinimumValue=15,
    MaximumValue=15,
    TimeTolerance=0.10000000000000001,
    LogBoundary="Left",
)
GroupWorkspaces(InputWorkspaces="30807_entry-Off_Off", OutputWorkspace="wsg")
GroupWorkspaces(InputWorkspaces="30807_entry-Off_Off", OutputWorkspace="wsg")
MagnetismReflectometryReduction(
    InputWorkspace="wsg",
    NormalizationWorkspace="30794_entry-Off_Off",
    SignalPeakPixelRange="181,195",
    SignalBackgroundPixelRange="49,88",
    NormPeakPixelRange="202,216",
    NormBackgroundPixelRange="94,104",
    LowResDataAxisPixelRange="69,172",
    LowResNormAxisPixelRange="71,175",
    TimeAxisRange="11413.6,45388.8",
    RoundUpPixel=False,
    SpecularPixel=187.80000000000001,
    FinalRebin=False,
    QMin=0.001,
    QStep=-0.02,
    TimeAxisStep=400,
    ConstQTrim=0.10000000000000001,
    OutputWorkspace="r30807",
)
Scale(InputWorkspace="r30807", OutputWorkspace="r30807", Factor=0.57943431658725553)
Scale(InputWorkspace="r30807", OutputWorkspace="r30807_scaled", Factor=0.61799999999999999)
AddSampleLog(Workspace="r30807_scaled", LogName="scaling_factor", LogText="0.618", LogType="Number")
# Run:30808
######################################################################
# Python Script Generated by GeneratePythonScript Algorithm
######################################################################
LoadEventNexus(Filename="/SNS/REF_M/IPTS-21391/nexus/REF_M_30796.nxs.h5", OutputWorkspace="raw_events")
FilterByLogValue(
    InputWorkspace="raw_events",
    OutputWorkspace="30796_entry-Off_Off",
    LogName="BL4A:SF:ICP:getDI",
    MinimumValue=15,
    MaximumValue=15,
    TimeTolerance=0.10000000000000001,
    LogBoundary="Left",
)
GroupWorkspaces(InputWorkspaces="30796_entry-Off_Off", OutputWorkspace="wsg")
LoadEventNexus(Filename="/SNS/REF_M/IPTS-21391/nexus/REF_M_30808.nxs.h5", OutputWorkspace="raw_events")
FilterByLogValue(
    InputWorkspace="raw_events",
    OutputWorkspace="30808_entry-Off_Off",
    LogName="BL4A:SF:ICP:getDI",
    MinimumValue=15,
    MaximumValue=15,
    TimeTolerance=0.10000000000000001,
    LogBoundary="Left",
)
GroupWorkspaces(InputWorkspaces="30808_entry-Off_Off", OutputWorkspace="wsg")
GroupWorkspaces(InputWorkspaces="30808_entry-Off_Off", OutputWorkspace="wsg")
GroupWorkspaces(InputWorkspaces="30808_entry-Off_Off", OutputWorkspace="wsg")
GroupWorkspaces(InputWorkspaces="30808_entry-Off_Off", OutputWorkspace="wsg")
GroupWorkspaces(InputWorkspaces="30808_entry-Off_Off", OutputWorkspace="wsg")
MagnetismReflectometryReduction(
    InputWorkspace="wsg",
    NormalizationWorkspace="30796_entry-Off_Off",
    SignalPeakPixelRange="183,197",
    SignalBackgroundPixelRange="49,88",
    NormPeakPixelRange="202,216",
    NormBackgroundPixelRange="94,104",
    LowResDataAxisPixelRange="69,172",
    LowResNormAxisPixelRange="82,175",
    TimeAxisRange="11413.6,45388.8",
    RoundUpPixel=False,
    SpecularPixel=189,
    FinalRebin=False,
    QMin=0.001,
    QStep=-0.02,
    TimeAxisStep=400,
    ConstQTrim=0.10000000000000001,
    OutputWorkspace="r30808",
)
Scale(InputWorkspace="r30808", OutputWorkspace="r30808", Factor=0.19100143649212772)
Scale(InputWorkspace="r30808", OutputWorkspace="r30808_scaled", Factor=0.72276980360217025)
AddSampleLog(Workspace="r30808_scaled", LogName="scaling_factor", LogText="0.722769803602", LogType="Number")
