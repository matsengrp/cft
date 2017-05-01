/*eslint-env browser*/
import React from "react";
import { connect } from "react-redux";
//import { loadJSONs } from (auspiceRoot + "actions/loadData");
import { loadJSONs } from "./auspice/src/actions/loadData";
import { RESET_CONTROLS, NEW_DATASET } from  "./auspice/src/actions/types";
//import { restoreStateFromURL, turnURLtoDataPath } from "./auspice/src/util/urlHelpers"
import "whatwg-fetch"; // setup polyfill
import Title from "./auspice/src/components/framework/title";
import Background from "./auspice/src/components/framework/background";
import ToggleSidebarTab from "./auspice/src/components/framework/toggle-sidebar-tab";
import Controls from "./auspice/src/components/controls/controls";
import Frequencies from "./auspice/src/components/charts/frequencies";
import Entropy from "./auspice/src/components/charts/entropy";
//import Map from "./auspice/src/components/map/map";
import TreeView from "./auspice/src/components/tree/treeView";
import queryString from "query-string";
import * as globals from "./auspice/src/util/globals";
import Sidebar from "react-sidebar";
import Flex from "./auspice/src/component/framework/flex";
import { titleStyles } from "./auspice/src/globalStyles";
import TitleBar from "./auspice/src/component/framework/title-bar";
//import Footer from "./auspice/src/component/framework/footer";
//import { analyticsNewPage } from "./auspice/src/util/googleAnalytics";



/* BRIEF REMINDER OF PROPS AVAILABLE TO APP:
  React-Router v4 injects length, action, location, push etc into props,
    but perhaps it's more consistent if we access these through
    this.context.router.
  Regardless, changes in URL will trigger the lifecycle methods
    here as that is a prop of this component, whether we use it or not
  see https://reacttraining.com/react-router
*/
@connect()
class App extends React.Component {
  constructor(props) {
    super(props);
    /* window listener to see when width changes cross thrhershold to toggle sidebar */
    /* A note on sidebar terminology:
    sidebarOpen (AFAIK) is only used via touch drag events
    sidebarDocked is the prop used on desktop.
    While these states could be moved to redux, they would need
    to be connected to here, triggering an app render anyways
    */
    const mql = window.matchMedia(`(min-width: ${globals.controlsHiddenWidth}px)`);
    mql.addListener(() => this.setState({sidebarDocked: this.state.mql.matches}));
    this.state = {
      mql,
      sidebarDocked: mql.matches,
      sidebarOpen: false
    };
    analyticsNewPage();
  }
  static propTypes = {
    dispatch: React.PropTypes.func.isRequired
  }
  static contextTypes = {
    router: React.PropTypes.object.isRequired
  }

  // componentWillMount() {
  // }

  componentWillMount() {
    /* parse URL, set URL, load data etc
    */

    this.props.dispatch({type: RESET_CONTROLS});
    // TODO Define our own versions of these functions
    //const data_path = turnURLtoDataPath(this.context.router);
    //restoreStateFromURL(this.context.router, this.props.dispatch);
    if (data_path) {
      this.props.dispatch({type: NEW_DATASET, data: this.context.router.history.location.pathname});
      this.props.dispatch(loadJSONs(data_path));
    } else {
      console.log("<app> couldn't work out the dataset to load. Bad.");
    }
  }

  componentDidUpdate() {
    /* back/forward buttons used (i.e. app doesn't reload, things don't
    remount, but this is the place to detect URL changes)
    */
    // // debugging hook?
    // console.log("app.js CDU")
    // console.log("redux datasetPathName:", this.props.datasetPathName);
    // this.maybeFetchDataset();
  }

  render() {
    return (
      <Sidebar
        sidebar={
          <div>
            <TitleBar minified={true}/>
            <Controls/>
            <ToggleSidebarTab
              open={this.state.sidebarDocked}
              handler={() => {this.setState({sidebarDocked: !this.state.sidebarDocked});}}
            />
          </div>
        }
        open={this.state.sidebarOpen}
        docked={this.state.sidebarDocked}
        onSetOpen={(a) => {this.setState({sidebarOpen: a});}}>
        <Background>
          {this.state.sidebarOpen || this.state.sidebarDocked ? <div/> :
            <ToggleSidebarTab
              open={this.state.sidebarDocked}
              handler={() => {this.setState({sidebarDocked: !this.state.sidebarDocked});}}
            />
          }
          <TreeView
            query={queryString.parse(this.context.router.history.location.search)}
            sidebar={this.state.sidebarOpen || this.state.sidebarDocked}
          />
          <Frequencies/>
          <Entropy
            sidebar={this.state.sidebarOpen || this.state.sidebarDocked}
          />
        </Background>
      </Sidebar>
    );
  }
}

export default App;

