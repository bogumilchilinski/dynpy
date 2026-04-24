from abc import ABC, abstractmethod
from typing import List, Optional, Any
from datetime import datetime
import requests
from typing import Any, Dict, Optional


class Project(ABC):
    """
    Abstract base class representing a generic issue-tracking project
    (GitHub, Redmine, Jira, etc.).

    This class defines a common interface for project backends.
    """
    
    _project_url = "https://pmt.example.com"
    _te_endpoint = "/time_entries.json"
    _api_key = None

    def __init__(self, project_name: Optional[str] = None,api_key: Optional[str] = None):
        self._project_name = project_name
        self._api_key = api_key
        self.__class__._api_key = api_key

    # ------------------------------------------------------------------
    # CONNECTION LIFECYCLE
    # ------------------------------------------------------------------

    #@abstractmethod
    def open(self) -> None:
        """Open connection to the backend system."""
        pass

    #@abstractmethod
    def close(self) -> None:
        """Close connection to the backend system."""
        pass

    # ------------------------------------------------------------------
    # PROJECT ACCESS
    # ------------------------------------------------------------------

    #@abstractmethod
    def get_project(self) -> Any:
        """Return backend-specific project object."""
        pass

    # ------------------------------------------------------------------
    # ISSUES
    # ------------------------------------------------------------------

    #@abstractmethod
    def get_issues(
        self,
        state: str = "all",
        assignee: Optional[str] = None,
        since: Optional[datetime] = None,
        sort: Optional[str] = None,
    ) -> List[Any]:
        """Returns a list of issues."""
        pass

    def get_time_entries(
        self,
        issue: Optional[str] = None,
        assignee: Optional[str] = None,
        since: Optional[datetime] = None,
        sort: Optional[str] = None,
    ) -> List[Any]:

    """Returns a list of time entries."""
    
    
    base_url = f"{project_url.rstrip('/')}{self._endpoint}"
    
    # parametry zapytania GET
    params = urlencode({"issue_id": issue})
    url = f"{base_url}?{params}"

    req_result = RedmineProject('',api_key=api_key)._request_data('GET',url)
    
    return req_result



    #@abstractmethod
    def get_issue(self, issue_id: int) -> Any:
        """Return a single issue by its ID or number."""
        pass

    # ------------------------------------------------------------------
    # ISSUE MODIFICATION
    # ------------------------------------------------------------------

    #@abstractmethod
    def create_issue(
        self,
        title: str,
        description: Optional[str] = None,
        assignee: Optional[str] = None,
        labels: Optional[List[str]] = None,
    ) -> Any:
        """Create a new issue."""
        pass

    #@abstractmethod
    def close_issue(self, issue_id: int) -> None:
        """Close an existing issue."""
        pass


    def _request_data(self,
        method: str,
        url: str,
        payload: Optional[Dict[str, Any]] = None,
        headers: Optional[Dict[str, Any]] = None,
        timeout: int = 30,
        ):
        
        response = requests.request(
            method=method.upper(),
            url=url,
            headers=headers,
            json=payload,
            timeout=timeout,
        )

        response.raise_for_status()

        if response.content:
            return response.json()
        
        return response



    def system_request(self,
        method: str,
        endpoint: str,
        payload: Optional[Dict[str, Any]] = None,
        system_url: str = "https://pmt.example.com",
        api_key: str = "",
        timeout: int = 30,
    ) -> Dict[str, Any]:
        """
        Execute a generic request to PMT REST API.

        Parameters
        ----------
        method : str
            HTTP method: 'GET', 'POST', 'PUT', 'DELETE'
        endpoint : str
            API endpoint, e.g. '/issues.json'
        payload : dict, optional
            JSON payload sent in request body
        system_url : str
            Base URL
        api_key : str
            
            PMT API key
        timeout : int
            Request timeout in seconds

        Returns
        -------
        dict
            Parsed JSON response (if any)
        """

        url = f"{system_url.rstrip('/')}/{endpoint.lstrip('/')}"
        headers = {
            "Content-Type": "application/json",
            #"api-key": api_key,
        }

        response = self._request_data(
            method=method.upper(),
            url=url,
            headers=headers,
            json=payload,
            timeout=timeout,
        )


        return response



class CodeRepository(Project):
    pass


class IssueEntry(ABC):
    """
    Abstract base class representing a generic entry within a Redmine issue.
    This serves as a common type for polymorphic behavior.
    """
    pass

class Comment(IssueEntry):
    """
    Represents a textual comment added to an issue.
    """
    def __init__(self, content: str, author: str):
        self.content = content
        self.author = author

class Note(IssueEntry):
    """
    Represents a journal note, which may include internal observations or updates.
    """
    def __init__(self, text: str, is_private: bool = False):
        self.text = text
        self.is_private = is_private

class TimeEntry(IssueEntry):
    """
    Represents a log of time spent working on a specific task.
    """
    def __init__(self, hours: float, activity: str):
        self.hours = hours
        self.activity = activity
