

from __future__ import annotations
from dataclasses import dataclass
from typing import Protocol, Optional, Sequence, Mapping, Any
import os, textwrap, json









class Summarizer(Protocol):
    def summarize(self, *, prompt: str, max_tokens: int = 700) -> str: ...


class HeuristicSummarizer:
    """Bardzo prosty fallback — wyciąga najdłuższe zdania + metryki."""
    def summarize(self, *, prompt: str, max_tokens: int = 700) -> str:
        sentences = [s.strip() for s in prompt.replace("\n", " ").split(".") if s.strip()]
        key_words = ("RMSE","MAE","MSE","ACC","F1","%","±","W","Wh","SOC","błąd","dokładność")
        scored = sorted(
            sentences,
            key=lambda s: sum(int(k in s) for k in key_words) + min(len(s),160)/160.0,
            reverse=True
        )
        return ". ".join(scored[:5]) + "."







@dataclass
class OpenAIChatSummarizer:
    model: str = "gpt-4o-mini"
    temperature: float = 0.2

    def _client(self):
        try:
            from openai import OpenAI  # pip install openai
        except Exception as e:
            raise RuntimeError("Brak pakietu openai. Zainstaluj: pip install openai") from e
        if not os.environ.get("OPENAI_API_KEY"):
            raise RuntimeError("Brak OPENAI_API_KEY w środowisku.")
        base_url = os.environ.get("OPENAI_BASE_URL") or None
        return OpenAI(base_url=base_url)

    def summarize(self, *, prompt: str, max_tokens: int = 700) -> str:
        client = self._client()
        resp = client.chat.completions.create(
            model=self.model,
            temperature=self.temperature,
            messages=[
                {"role": "system", "content": "Jesteś asystentem akademickim. Podsumowujesz wyniki eksperymentów po polsku, z naciskiem na metryki i wnioski."},
                {"role": "user", "content": prompt},
            ],
            max_tokens=max_tokens,
        )
        return (resp.choices[0].message.content or "").strip()








@dataclass
class OllamaSummarizer:

    model: str = "llama3"
    host: str = "http://127.0.0.1:11434"
    temperature: float = 0.2
    top_p: float = 0.95

    def _build_prompt(self, prompt: str) -> str:
        system = ("Jesteś asystentem akademickim. Podsumowujesz wyniki eksperymentów po polsku, "
                  "z naciskiem na liczby, metryki i wnioski.")
        return f"[SYSTEM]\n{system}\n\n[USER]\n{prompt}\n\n[ASSISTANT]"

    def summarize(self, *, prompt: str, max_tokens: int = 700) -> str:
        import urllib.request, urllib.error
        url = f"{self.host}/api/generate"
        payload = {
            "model": self.model,
            "prompt": self._build_prompt(prompt),
            "options": {
                "temperature": self.temperature,
                "top_p": self.top_p,
                "num_predict": max(64, min(300, max_tokens)),
            },
            "stream": False,
        }
        data = json.dumps(payload).encode("utf-8")
        req = urllib.request.Request(url, data=data, headers={"Content-Type": "application/json"}, method="POST")
        with urllib.request.urlopen(req, timeout=120) as resp:
            raw = resp.read().decode("utf-8")
            obj = json.loads(raw) if raw else {}
        out = obj.get("response") if isinstance(obj, dict) else None
        if isinstance(out, str) and out.strip():
            return out.strip()
        raise RuntimeError(f"Ollama: nieoczekiwana odpowiedź: {obj}")









def _format_metrics(metrics: Optional[Mapping[str, Any]]) -> str:
    if not metrics:
        return "—"
    try:
        return "\n".join(f"- {k}: {metrics[k]}" for k in sorted(metrics))
    except Exception:
        return json.dumps(metrics, ensure_ascii=False, indent=2)

def build_summary_prompt(*, objective: str, experimental_setup: str,
                         results_table=None, metrics: Optional[Mapping[str, Any]]=None,
                         observations: Optional[Sequence[str]]=None,
                         constraints: str="maks. 130–400 słów, język polski, styl raportowy") -> str:
    obs_text = "\n".join(f"- {o}" for o in (observations or []))
    table_txt = "brak tabeli wyników"
    if results_table is not None:
        try:
            head = results_table.head(15).to_string(index=False)
            table_txt = f"(pierwsze wiersze)\n{head}"
        except Exception:
            pass
    return textwrap.dedent(f"""
    Cel eksperymentu:
    {objective}

    Opis stanowiska / warunków:
    {experimental_setup}

    Streszczenie tabeli wyników:
    {table_txt}

    Kluczowe metryki:
    {_format_metrics(metrics)}

    Dodatkowe obserwacje:
    {obs_text or "—"}

    Zadanie: Sporządź zwarte podsumowanie wniosków z naciskiem na liczby, trend i ograniczenia. {constraints}.
    """).strip()

def summarize_experiment(*, summarizer: Summarizer, objective: str, experimental_setup: str,
                         results_table=None, metrics: Optional[Mapping[str, Any]]=None,
                         observations: Optional[Sequence[str]]=None, max_tokens: int=700) -> str:
    prompt = build_summary_prompt(
        objective=objective,
        experimental_setup=experimental_setup,
        results_table=results_table,
        metrics=metrics,
        observations=observations,
    )
    return summarizer.summarize(prompt=prompt, max_tokens=max_tokens)





__all__ = [
    "Summarizer",
    "HeuristicSummarizer",
    "OpenAIChatSummarizer",
    "OllamaSummarizer",
    "build_summary_prompt",
    "summarize_experiment",
    "report_text_with_citations",
]
