from time import time

def myTimer(state,data={}):

    # add elapsed time to the current state of the timer
    if "currentState" in data:
        currentState = data["currentState"]
        data[currentState] += time() - data["lastTime"]

    if state == "summary":
        outtext = "Summary of time taken:\n"
        total = 0
        for state in data["allStates"]:
            if not state == "lastTime" and not state == "currentState":
                total += data[state]
        for state in data["allStates"]:
            if not state == "lastTime" and not state == "currentState":
                outtext += "%-31s: %10.3f seconds, %10.3f minutes, %10.2f percent\n" % (state,data[state],data[state]/60,100*data[state]/total)
        outtext += "%-31s: %10.3f seconds, %10.3f minutes\n" % ("Total",total,total/60)
        return outtext
    elif not state in data:
        data[state] = 0
        # keep track of states and the order in which they were seen
        if "allStates" in data:
            data["allStates"].append(state)
        else:
            data["allStates"] = [state]

    # change to the state just starting now
    data["currentState"] = state
    data["lastTime"] = time()

    return data
