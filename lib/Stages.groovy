import Constants

class Stages {

  public static set_stages(mode_str, log) {
    def mode_enum = Constants.WorkflowType.valueOf(mode_str.toUpperCase())
    def stages = []
    switch(mode_enum) {
      case Constants.WorkflowType.FULL:
        stages = Constants.Stage.values()
        break
      case Constants.WorkflowType.GRIDSS_PURPLE_LINX:
        stages = [
          Constants.Stage.AMBER,
          Constants.Stage.COBALT,
          Constants.Stage.GRIDSS,
          Constants.Stage.GRIPSS,
          Constants.Stage.LINX,
          Constants.Stage.PURPLE,
          Constants.Stage.SVPREP,
        ]
        break
      case Constants.WorkflowType.GRIDSS:
        stages = [
          Constants.Stage.GRIDSS,
          Constants.Stage.GRIPSS,
          Constants.Stage.SVPREP,
        ]
        break
      case Constants.WorkflowType.PURPLE:
        stages = [
          Constants.Stage.PURPLE,
        ]
        break
      case Constants.WorkflowType.LINX:
        stages = [
          Constants.Stage.LINX,
        ]
        break
      case Constants.WorkflowType.LILAC:
        stages = [
          Constants.Stage.LILAC
        ]
        break
      case Constants.WorkflowType.TEAL:
        stages = [
          Constants.Stage.TEAL
        ]
        break
      default:
        log.error "\nERROR: recieved invalid mode '${mode_str}'"
        System.exit(1)
    }
    return stages
  }

}
