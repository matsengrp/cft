import React from "react";
import ReactDOM from "react-dom";
import ReactGA from "react-ga";
import { Provider } from "react-redux";
import { BrowserRouter, Route, Switch } from "react-router-dom";
import injectTapEventPlugin from "react-tap-event-plugin";
import configureStore from "./auspice/src/store";
import App from "./app";
import BrowserDimensionMonitor from "./auspice/src/components/framework/browserDimensionMonitor";
import "./auspice/src/css/global.css";
import "./auspice/src/css/browserCompatability.css";
import "./auspice/src/css/bootstrapCustomized.css";
import "./auspice/src/css/datePicker.css";
import "./auspice/src/css/static.css";
//import { outboundLinkWithAnalytics } from "./auspice/src/util/googleAnalytics";


class Root extends React.Component {
  render() {
    return (
      <Provider store={store}>
        <BrowserRouter>
          <div>
            <BrowserDimensionMonitor/>
            <Switch>
              <Route path="/dengue*" component={App}/>
              <Route path="*" component={App}/>
            </Switch>
          </div>
        </BrowserRouter>
      </Provider>
    );
  }
}

/*  to fix iOS's dreaded 300ms tap delay, we need this plugin
NOTE Facebook is not planning on supporting tap events (#436)
because browsers are fixing/removing the click delay.
Unfortunately it will take a lot of time before all mobile
browsers (including iOS' UIWebView) will and can be updated.
https://github.com/zilverline/react-tap-event-plugin
*/
injectTapEventPlugin();

ReactDOM.render(<Root/>, document.getElementById("react-root"));

